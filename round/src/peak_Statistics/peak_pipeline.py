"""Peak-only analysis of P(phi,psi | m-1, m, m+1).

For every central residue b (20 groups):
  1. extract peaks of each of the 400 contexts (smoothed, sigma=1.5)
  2. cluster peak positions across contexts -> shared peak dictionary
  3. compute per-context peak-mass weights over the shared peaks
  4. compare the marginal P(phi,psi|b) peaks with the dictionary
     (matched / swallowed / new)
  5. flag outlier contexts whose peaks cannot be explained by the dictionary
     and copy their images from result/static/img_all

Outputs: result/peak/{peaks,weights,marginal}/..., outliers.csv, report.txt,
         images/<center>/*.png (copied, never re-rendered).
"""
import json
import os
import shutil

import numpy as np

NB = 72
BW = 5.0
ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
ROUND = r"D:\CodingPrograms\Protein Design\round"
PROBS = os.path.join(ROUND, "result", "static", "triplet_probs.npy")
IMG_ALL = os.path.join(ROUND, "result", "static", "img_all")
OUT = os.path.join(ROUND, "result", "peak")

SIGMA = 1.5
PEAK_FRAC = 0.08      # peak must reach 5% of the context maximum
MAX_PEAKS = 6
CLUSTER_R = 3.0       # bins (10 deg) on the torus
MIN_SUPPORT = 10       # contexts
MIN_MASS_FRAC = 0.02  # cluster must carry >=1% of the total peak mass
OUTLIER_UNMATCHED = 0.40
MIN_OCC = 100


def ksig(n, s):
    i = np.arange(n, dtype=float)
    d = np.minimum(i, n - i)
    k = np.exp(-0.5 * (d / s) ** 2)
    return k / k.sum()


def smooth(m, K):
    out = np.fft.ifft2(np.fft.fft2(m) * K).real
    np.clip(out, 0.0, None, out=out)
    s = out.sum()
    return out / s if s > 0 else out


def local_maxima(m):
    lm = np.ones_like(m, dtype=bool)
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            if dx == 0 and dy == 0:
                continue
            lm &= m > np.roll(np.roll(m, dx, axis=0), dy, axis=1)
    return lm


def torus_dist(i1, j1, i2, j2):
    di = abs(i1 - i2)
    dj = abs(j1 - j2)
    di = min(di, NB - di)
    dj = min(dj, NB - dj)
    return float(np.hypot(di, dj))


def window_mass(m, i, j, r=2):
    total = 0.0
    for di in range(-r, r + 1):
        for dj in range(-r, r + 1):
            if di * di + dj * dj <= r * r:
                total += m[(i + di) % NB, (j + dj) % NB]
    return float(total)


def extract_peaks(m):
    lm = local_maxima(m)
    thr = PEAK_FRAC * m.max()
    cand = np.argwhere(lm & (m >= thr))
    if len(cand) == 0:
        return []
    vals = m[cand[:, 0], cand[:, 1]]
    order = np.argsort(-vals)[:MAX_PEAKS]
    peaks = []
    for k in order:
        i, j = int(cand[k, 0]), int(cand[k, 1])
        peaks.append({
            "i": i,
            "j": j,
            "phi": -180.0 + BW * i + BW / 2,
            "psi": -180.0 + BW * j + BW / 2,
            "height": float(m[i, j]),
            "mass": window_mass(m, i, j, 2),
        })
    return peaks


def main():
    os.makedirs(os.path.join(OUT, "peaks"), exist_ok=True)
    os.makedirs(os.path.join(OUT, "weights"), exist_ok=True)
    os.makedirs(os.path.join(OUT, "marginal"), exist_ok=True)
    os.makedirs(os.path.join(OUT, "images"), exist_ok=True)

    probs = np.load(PROBS, mmap_mode="r")
    R = np.asarray(probs, dtype=np.float64).reshape(20, 20, 20, NB, NB)
    occ = (R > 0).sum(axis=(3, 4)).astype(float)

    k = ksig(NB, SIGMA)
    K = np.fft.fft(k)[:, None] * np.fft.fft(k)[None, :]
    S = np.empty_like(R)
    for a in range(20):
        for b in range(20):
            for c in range(20):
                S[a, b, c] = smooth(R[a, b, c], K)

    report = []
    report.append("peak-only analysis")
    report.append(f"sigma={SIGMA} peak>={PEAK_FRAC:.0%} max, max {MAX_PEAKS} peaks/context, "
                  f"cluster radius {CLUSTER_R} bins, min support {MIN_SUPPORT}")

    outliers = []
    summary = []

    for b in range(20):
        # --- collect peaks of the 400 contexts -----------------------------
        all_peaks = []
        for a in range(20):
            for c in range(20):
                w = occ[a, b, c]
                if w < MIN_OCC:
                    continue
                for p in extract_peaks(S[a, b, c]):
                    if p["mass"] <= 0:
                        continue
                    p.update({"a": a, "c": c, "occ": w})
                    all_peaks.append(p)

        if not all_peaks:
            summary.append({"center": ALPHABET[b], "clusters": 0})
            continue

        # --- greedy clustering on the torus --------------------------------
        all_peaks.sort(key=lambda p: -p["mass"])
        clusters = []
        for p in all_peaks:
            best = None
            best_d = CLUSTER_R + 1e-9
            for cl in clusters:
                d = torus_dist(p["i"], p["j"], cl["i"], cl["j"])
                if d < best_d:
                    best_d = d
                    best = cl
            if best is None:
                clusters.append({"i": p["i"], "j": p["j"], "mass": p["mass"],
                                 "support": 1, "members": [p], "sum_w": p["mass"]})
            else:
                best["mass"] += p["mass"]
                best["support"] += 1
                best["members"].append(p)

        total_mass = sum(cl["mass"] for cl in clusters)
        clusters.sort(key=lambda cl: -cl["mass"])
        kept = []
        cum = 0.0
        for cl in clusters:
            frac = cl["mass"] / total_mass
            if cl["support"] < MIN_SUPPORT or frac < MIN_MASS_FRAC:
                continue
            # weighted circular-ish mean (positions are bin indices, no wrap here)
            wsum = sum(m["mass"] for m in cl["members"])
            ci = sum(m["i"] * m["mass"] for m in cl["members"]) / wsum
            cj = sum(m["j"] * m["mass"] for m in cl["members"]) / wsum
            sd_i = np.sqrt(sum(m["mass"] * (m["i"] - ci) ** 2 for m in cl["members"]) / wsum)
            sd_j = np.sqrt(sum(m["mass"] * (m["j"] - cj) ** 2 for m in cl["members"]) / wsum)
            cl.update({"ci": ci, "cj": cj, "sd_i": float(sd_i), "sd_j": float(sd_j),
                       "frac": frac})
            kept.append(cl)
            cum += frac
            if cum >= 0.90:
                break

        # --- per-context weights / unmatched mass --------------------------
        weight = np.zeros((20, 20, len(kept)))
        ctx_matched = []
        for a in range(20):
            for c in range(20):
                if occ[a, b, c] < MIN_OCC:
                    continue
                ctx_peaks = extract_peaks(S[a, b, c])
                tm = sum(p["mass"] for p in ctx_peaks) or 1.0
                matched = 0.0
                for p in ctx_peaks:
                    best_k = -1
                    best_d = CLUSTER_R + 1e-9
                    for kk, cl in enumerate(kept):
                        d = torus_dist(p["i"], p["j"], cl["ci"], cl["cj"])
                        if d < best_d:
                            best_d = d
                            best_k = kk
                    if best_k >= 0:
                        weight[a, c, best_k] += p["mass"] / tm
                        matched += p["mass"] / tm
                unmatched = 1.0 - matched
                ctx_matched.append(matched)
                if unmatched > OUTLIER_UNMATCHED:
                    outliers.append({
                        "pdb": "", "a": ALPHABET[a], "b": ALPHABET[b], "c": ALPHABET[c],
                        "unmatched": round(float(unmatched), 3),
                        "occ": int(occ[a, b, c]),
                    })
        weight = np.clip(weight, 0.0, None)
        np.save(os.path.join(OUT, "weights", f"center_{ALPHABET[b]}.npy"), weight)

        # --- marginal peaks ------------------------------------------------
        marg = np.zeros((NB, NB))
        wsum = 0.0
        for a in range(20):
            for c in range(20):
                if occ[a, b, c] >= MIN_OCC:
                    marg += occ[a, b, c] * S[a, b, c]
                    wsum += occ[a, b, c]
        marg = marg / wsum if wsum > 0 else marg
        np.save(os.path.join(OUT, "marginal", f"center_{ALPHABET[b]}.npy"), marg)
        mpeaks = extract_peaks(marg)

        # swallowed: shared-peak mass that has no marginal peak nearby
        matched_cluster_mass = 0.0
        for cl in kept:
            if any(torus_dist(p["i"], p["j"], cl["ci"], cl["cj"]) <= CLUSTER_R
                   for p in mpeaks):
                matched_cluster_mass += cl["frac"]
        swallowed = 1.0 - matched_cluster_mass

        entry = {
            "center": ALPHABET[b],
            "n_clusters": len(kept),
            "context_peak_coverage": round(float(np.mean(ctx_matched)), 3) if ctx_matched else 0.0,
            "clusters": [{
                "phi": round(-180.0 + BW * cl["ci"] + BW / 2, 2),
                "psi": round(-180.0 + BW * cl["cj"] + BW / 2, 2),
                "mass_frac": round(float(cl["frac"]), 3),
                "support": int(cl["support"]),
                "sd_phi_deg": round(float(cl["sd_i"]) * BW, 1),
                "sd_psi_deg": round(float(cl["sd_j"]) * BW, 1),
            } for cl in kept],
            "marginal_peaks": [{
                "phi": round(p["phi"], 1), "psi": round(p["psi"], 1),
                "height": round(p["height"], 5),
            } for p in mpeaks],
            "marginal_peak_mass_swallowed": round(float(swallowed), 3),
        }
        summary.append(entry)
        with open(os.path.join(OUT, "peaks", f"center_{ALPHABET[b]}.json"), "w",
                  encoding="utf-8") as f:
            json.dump(entry, f, ensure_ascii=False, indent=1)

    # --- outliers csv + copied images -------------------------------------
    seen = set()
    with open(os.path.join(OUT, "outliers.csv"), "w", encoding="utf-8") as f:
        f.write("a,b,c,unmatched,occ,image\n")
        for o in outliers:
            key = (o["a"], o["b"], o["c"])
            if key in seen:
                continue
            seen.add(key)
            src = os.path.join(IMG_ALL, o["b"], f"{o['a']}-{o['b']}-{o['c']}.png")
            dst_dir = os.path.join(OUT, "images", o["b"])
            os.makedirs(dst_dir, exist_ok=True)
            dst = os.path.join(dst_dir, f"{o['a']}-{o['b']}-{o['c']}.png")
            copied = ""
            if os.path.exists(src):
                shutil.copy2(src, dst)
                copied = os.path.relpath(dst, OUT)
            f.write(f"{o['a']},{o['b']},{o['c']},{o['unmatched']},{o['occ']},{copied}\n")

    # --- report ------------------------------------------------------------
    report.append("")
    report.append(f"{'ctr':>3} {'K':>3} {'cover':>6} {'swallow':>8} {'clusters (phi,psi,mass,support)'}")
    for e in summary:
        if "clusters" not in e:
            report.append(f"{e['center']:>3}   0")
            continue
        cl = "  ".join(
            f"({c['phi']:.0f},{c['psi']:.0f},{c['mass_frac']:.2f},{c['support']})"
            for c in e["clusters"][:5])
        report.append(f"{e['center']:>3} {e['n_clusters']:>3} {e['context_peak_coverage']:>6.2f} "
                      f"{e['marginal_peak_mass_swallowed']:>8.3f}  {cl}")
    report.append("")
    report.append(f"outlier contexts: {len(seen)}")
    text = "\n".join(report)
    print(text)
    with open(os.path.join(OUT, "report.txt"), "w", encoding="utf-8") as f:
        f.write(text + "\n")


if __name__ == "__main__":
    main()

