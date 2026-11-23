"""Export the final peak representation into result/peak/fin_peaks.

Structure:
  dictionary.json   shared peak dictionary per center residue
                    { "A": [ {id, phi, psi, mass_frac, support, sd_phi, sd_psi}, ... ] }
  contexts.jsonl    one line per context (a,b,c):
                    { a, b, c, matched: {peak_id: weight, ...},
                      extra: [ {phi, psi, weight}, ... ], coverage }
                    - matched peaks reference the dictionary by id (positions fixed)
                    - extra peaks carry their own position (the "outlier" peaks)
  outliers.jsonl    subset of contexts that have at least one extra peak,
                    with the extra peaks and their weights
  README.md
"""
import json
import os

import numpy as np

NB = 72
BW = 5.0
ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
ROUND = r"D:\CodingPrograms\Protein Design\round"
PROBS = os.path.join(ROUND, "result", "static", "triplet_probs.npy")
OUT = os.path.join(ROUND, "result", "peak", "fin_peaks")

SIGMA = 1.5
PEAK_FRAC = 0.08
MAX_PEAKS = 6
CLUSTER_R = 3.0
MIN_SUPPORT = 10
MIN_MASS_FRAC = 0.02
MIN_OCC = 100
EXTRA_MIN_W = 0.02   # only keep extras above 2% of the context peak mass


def ksig(n, s):
    i = np.arange(n, dtype=float)
    d = np.minimum(i, n - i)
    k = np.exp(-0.5 * (d / s) ** 2)
    return k / k.sum()


def local_maxima(m):
    lm = np.ones_like(m, dtype=bool)
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            if dx == 0 and dy == 0:
                continue
            lm &= m > np.roll(np.roll(m, dx, 0), dy, 1)
    return lm


def torus_dist(i1, j1, i2, j2):
    di = min(abs(i1 - i2), NB - abs(i1 - i2))
    dj = min(abs(j1 - j2), NB - abs(j1 - j2))
    return float(np.hypot(di, dj))


def extract(m):
    lm = local_maxima(m)
    thr = PEAK_FRAC * m.max()
    cand = np.argwhere(lm & (m >= thr))
    if len(cand) == 0:
        # fallback for very flat/sparse contexts: use the global maximum
        i, j = np.unravel_index(int(np.argmax(m)), m.shape)
        cand = np.array([[i, j]])
    vals = m[cand[:, 0], cand[:, 1]]
    order = np.argsort(-vals)[:MAX_PEAKS]
    peaks = []
    for k in order:
        i, j = int(cand[k, 0]), int(cand[k, 1])
        mass = 0.0
        for di in range(-2, 3):
            for dj in range(-2, 3):
                if di * di + dj * dj <= 4:
                    mass += m[(i + di) % NB, (j + dj) % NB]
        if mass > 0:
            peaks.append({"i": i, "j": j, "mass": float(mass)})
    return peaks


def main():
    os.makedirs(OUT, exist_ok=True)
    probs = np.load(PROBS, mmap_mode="r")
    R = np.asarray(probs, dtype=np.float64).reshape(20, 20, 20, NB, NB)
    occ = (R > 0).sum(axis=(3, 4)).astype(float)

    k = ksig(NB, SIGMA)
    K = np.fft.fft(k)[:, None] * np.fft.fft(k)[None, :]
    S = np.empty_like(R)
    for a in range(20):
        for b in range(20):
            for c in range(20):
                m = np.fft.ifft2(np.fft.fft2(R[a, b, c]) * K).real
                m = np.clip(m, 0.0, None)
                S[a, b, c] = m / m.sum() if m.sum() > 0 else m

    dictionary = {}
    ctx_lines = []
    outlier_lines = []

    for b in range(20):
        all_peaks = []
        for a in range(20):
            for c in range(20):
                if occ[a, b, c] < MIN_OCC:
                    continue
                for p in extract(S[a, b, c]):
                    p.update({"a": a, "c": c})
                    all_peaks.append(p)
        all_peaks.sort(key=lambda p: -p["mass"])

        clusters = []
        for p in all_peaks:
            best, bd = None, CLUSTER_R + 1e-9
            for cl in clusters:
                d = torus_dist(p["i"], p["j"], cl["i"], cl["j"])
                if d < bd:
                    bd, best = d, cl
            if best is None:
                clusters.append({"i": p["i"], "j": p["j"], "mass": p["mass"],
                                 "support": 1, "members": [p]})
            else:
                best["mass"] += p["mass"]
                best["support"] += 1
                best["members"].append(p)

        total = sum(cl["mass"] for cl in clusters) or 1.0
        clusters.sort(key=lambda cl: -cl["mass"])
        kept, cum = [], 0.0
        for cl in clusters:
            frac = cl["mass"] / total
            if cl["support"] < MIN_SUPPORT or frac < MIN_MASS_FRAC:
                continue
            wsum = sum(m["mass"] for m in cl["members"])
            ci = sum(m["i"] * m["mass"] for m in cl["members"]) / wsum
            cj = sum(m["j"] * m["mass"] for m in cl["members"]) / wsum
            cl.update({"ci": ci, "cj": cj, "frac": frac})
            kept.append(cl)
            cum += frac
            if cum >= 0.90:
                break

        dictionary[ALPHABET[b]] = [{
            "id": f"{ALPHABET[b]}{idx + 1}",
            "phi": round(-180.0 + BW * cl["ci"] + BW / 2, 1),
            "psi": round(-180.0 + BW * cl["cj"] + BW / 2, 1),
            "mass_frac": round(float(cl["frac"]), 3),
            "support": int(cl["support"]),
        } for idx, cl in enumerate(kept)]

        for a in range(20):
            for c in range(20):
                peaks = extract(S[a, b, c])
                tm = sum(p["mass"] for p in peaks) or 1.0
                matched = {}
                extras = []
                for p in peaks:
                    best_k, bd = -1, CLUSTER_R + 1e-9
                    for kk, cl in enumerate(kept):
                        d = torus_dist(p["i"], p["j"], cl["ci"], cl["cj"])
                        if d < bd:
                            bd, best_k = d, kk
                    w = p["mass"] / tm
                    if best_k >= 0:
                        pid = f"{ALPHABET[b]}{best_k + 1}"
                        matched[pid] = round(matched.get(pid, 0.0) + w, 4)
                    elif w >= EXTRA_MIN_W:
                        extras.append({
                            "phi": round(-180.0 + BW * p["i"] + BW / 2, 1),
                            "psi": round(-180.0 + BW * p["j"] + BW / 2, 1),
                            "weight": round(w, 4),
                        })
                coverage = sum(matched.values())
                rec = {
                    "a": ALPHABET[a],
                    "b": ALPHABET[b],
                    "c": ALPHABET[c],
                    "occ": int(occ[a, b, c]),
                    "matched": matched,
                    "extra": extras,
                    "coverage": round(float(coverage), 3),
                }
                ctx_lines.append(json.dumps(rec, ensure_ascii=False))
                major_extra = any(e["weight"] >= 0.20 for e in extras)
                if major_extra or coverage < 0.75:
                    outlier_lines.append(json.dumps(rec, ensure_ascii=False))

    with open(os.path.join(OUT, "dictionary.json"), "w", encoding="utf-8") as f:
        json.dump(dictionary, f, ensure_ascii=False, indent=1)
    with open(os.path.join(OUT, "contexts.jsonl"), "w", encoding="utf-8") as f:
        f.write("\n".join(ctx_lines) + "\n")
    with open(os.path.join(OUT, "outliers.jsonl"), "w", encoding="utf-8") as f:
        f.write("\n".join(outlier_lines) + "\n")

    readme = """# fin_peaks —— 峰-位置-权重 + 异类列表

- `dictionary.json`：每个中心残基的共享峰字典，峰以 `id` 引用，位置固定。
- `contexts.jsonl`：8000 个上下文，一行一条：
  - `matched`：`{峰id: 权重}`，权重是该峰占本上下文峰质量的比重；
  - `extra`：字典覆盖不到的峰，自带位置 `(phi, psi)` 和权重；
  - `coverage`：`matched` 权重之和。
- `outliers.jsonl`：`contexts.jsonl` 中带 `extra` 的子集，即异类列表。

异类不再单独维护：它们就是 `extra` 非空的上下文，其余字典峰权重自然为 0。
"""
    with open(os.path.join(OUT, "README.md"), "w", encoding="utf-8") as f:
        f.write(readme)

    print(f"contexts={len(ctx_lines)} outliers={len(outlier_lines)}")
    print("dictionary sizes:", {k2: len(v) for k2, v in dictionary.items()})


if __name__ == "__main__":
    main()
