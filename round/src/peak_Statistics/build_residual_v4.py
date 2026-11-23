"""v4: add context-local peaks so significant empirical clusters always get fitted.

For each context:
  - start from dictionary peaks + existing extras
  - find local maxima of the smoothed empirical histogram
  - add a local peak when it carries >=1% of mass and is >2.5 bins from any
    existing peak (max 4 per context)
  - exclusive assignment -> actual peak masses, residual base
  - wash small residual specks (5x5 neighbourhood < 0.8%)
Outputs: ground/residual_v4, fin/weights_mass_v4, fin/base_mass_v4,
         fin/extras_mass_v4.jsonl
"""
import json
import os

import numpy as np

NB = 72
BW = 5.0
ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
ROUND = r"D:\CodingPrograms\Protein Design\round"
PROBS = os.path.join(ROUND, "result", "static", "triplet_probs.npy")
FIN = os.path.join(ROUND, "result", "fin")
GROUND = os.path.join(ROUND, "result", "ground")
ASSIGN_R = 3.0
SIGMA = 1.5
LOCAL_FRAC = 0.01
LOCAL_DIST = 2.5
MAX_LOCAL = 4
WASH = 0.008


def ksig(n, s):
    i = np.arange(n, dtype=float)
    d = np.minimum(i, n - i)
    k = np.exp(-0.5 * (d / s) ** 2)
    return k / k.sum()


def torus_d(i, j, pi, pj):
    di = min(abs(i - pi), NB - abs(i - pi))
    dj = min(abs(j - pj), NB - abs(j - pj))
    return np.hypot(di, dj)


def torus_grid(pi, pj):
    i = np.arange(NB)
    di = np.minimum(np.abs(i - pi), NB - np.abs(i - pi))
    dj = np.minimum(np.abs(i - pj), NB - np.abs(i - pj))
    return np.sqrt(di[:, None] ** 2 + dj[None, :] ** 2)


def main():
    os.makedirs(os.path.join(GROUND, "residual_v4"), exist_ok=True)
    os.makedirs(os.path.join(FIN, "weights_mass_v4"), exist_ok=True)
    os.makedirs(os.path.join(FIN, "base_mass_v4"), exist_ok=True)

    with open(os.path.join(FIN, "dictionary.json"), encoding="utf-8") as f:
        dictionary = json.load(f)
    extras0 = {}
    with open(os.path.join(FIN, "extras_mass.jsonl"), encoding="utf-8") as f:
        for line in f:
            rec = json.loads(line)
            extras0.setdefault((rec["a"], rec["b"], rec["c"]), []).append(rec)

    probs = np.load(PROBS, mmap_mode="r")
    k1 = ksig(NB, 1.0)
    K1 = np.fft.fft(k1)[:, None] * np.fft.fft(k1)[None, :]
    k2 = ksig(NB, SIGMA)
    K2 = np.fft.fft(k2)[:, None] * np.fft.fft(k2)[None, :]
    box = np.zeros((NB, NB))
    for di in range(-2, 3):
        for dj in range(-2, 3):
            box[di % NB, dj % NB] = 1.0
    KB = np.fft.fft2(box)

    out_lines = []
    for b in range(20):
        letter = ALPHABET[b]
        K = len(dictionary[letter])
        res = np.zeros((20, 20, NB, NB), dtype=np.float32)
        wmass = np.zeros((20, 20, K), dtype=np.float32)
        bmass = np.zeros((20, 20), dtype=np.float32)

        for a in range(20):
            for c in range(20):
                idx = a * 400 + b * 20 + c
                emp = np.asarray(probs[idx], dtype=np.float64)
                s = emp.sum()
                if s > 0:
                    emp /= s

                plist = [(dictionary[letter][q]["phi"], dictionary[letter][q]["psi"]) for q in range(K)]
                e0 = extras0.get((ALPHABET[a], letter, ALPHABET[c]), [])
                plist += [(e["phi"], e["psi"]) for e in e0]

                # local maxima of the smoothed empirical histogram
                F = np.fft.ifft2(np.fft.fft2(emp) * K1).real
                lm = np.ones((NB, NB), bool)
                for di in (-1, 0, 1):
                    for dj in (-1, 0, 1):
                        if di or dj:
                            lm &= F > np.roll(np.roll(F, di, 0), dj, 1)
                cand = [(int(i), int(j)) for i, j in zip(*np.where(lm)) if F[i, j] >= LOCAL_FRAC * F.max()]
                cand.sort(key=lambda t: -F[t[0], t[1]])
                local = []
                for i, j in cand:
                    if len(local) >= MAX_LOCAL:
                        break
                    far = all(torus_d(i, j, int((p[0] + 180) / BW) % NB, int((p[1] + 180) / BW) % NB) > LOCAL_DIST
                              for p in plist + local)
                    if far:
                        local.append((i, j))
                plist += [(-180 + BW * i + BW / 2, -180 + BW * j + BW / 2) for i, j in local]

                if not plist:
                    res[a, c] = emp.astype(np.float32)
                    bmass[a, c] = 1.0
                    continue

                dists = np.stack([torus_grid(int((p[0] + 180) / BW) % NB, int((p[1] + 180) / BW) % NB)
                                  for p in plist])
                nearest = np.argmin(dists, axis=0)
                assigned = np.min(dists, axis=0) <= ASSIGN_R

                for kk in range(K):
                    wmass[a, c, kk] = emp[(nearest == kk) & assigned].sum()
                for kk in range(K, len(plist)):
                    m = float(emp[(nearest == kk) & assigned].sum())
                    if kk >= K + len(e0):
                        phi, psi = plist[kk]
                    else:
                        phi, psi = e0[kk - K]["phi"], e0[kk - K]["psi"]
                    out_lines.append(json.dumps({"a": ALPHABET[a], "b": letter, "c": ALPHABET[c],
                                                 "phi": round(float(phi), 2), "psi": round(float(psi), 2),
                                                 "weight": round(m, 6)}, ensure_ascii=False))

                residual = np.where(assigned, 0.0, emp)
                bmass[a, c] = float(residual.sum())
                r = np.fft.ifft2(np.fft.fft2(residual) * K2).real
                r = np.clip(r, 0.0, None)
                nb = np.fft.ifft2(np.fft.fft2(r) * KB).real
                r = np.where(nb >= WASH, r, 0.0)
                r = r / r.sum() if r.sum() > 0 else np.ones((NB, NB)) / (NB * NB)
                res[a, c] = r.astype(np.float32)

        np.save(os.path.join(GROUND, "residual_v4", f"center_{letter}.npy"), res)
        np.save(os.path.join(FIN, "weights_mass_v4", f"center_{letter}.npy"), wmass)
        np.save(os.path.join(FIN, "base_mass_v4", f"center_{letter}.npy"), bmass)
        print("v4", letter)

    with open(os.path.join(FIN, "extras_mass_v4.jsonl"), "w", encoding="utf-8") as f:
        f.write("\n".join(out_lines) + "\n")
    print("extras rows", len(out_lines))


if __name__ == "__main__":
    main()
