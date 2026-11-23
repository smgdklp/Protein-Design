"""Rebuild weights and base with exclusive assignment.

For each context:
  - every cell is assigned to its nearest peak (dictionary or extra) within R bins
  - peak weights become ACTUAL probability mass (not fractions of peak mass)
  - residual (unassigned) mass is smoothed -> context-specific base
Outputs:
  result/ground/residual/center_X.npy   (20,20,72,72) normalized residual base
  result/fin/weights_mass/center_X.npy  (20,20,K)     actual peak mass
  result/fin/base_mass/center_X.npy     (20,20)       residual mass
  result/fin/extras_mass.jsonl          extra peaks with actual mass
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


def ksig(n, s):
    i = np.arange(n, dtype=float)
    d = np.minimum(i, n - i)
    k = np.exp(-0.5 * (d / s) ** 2)
    return k / k.sum()


def torus_dist_grid(pi, pj):
    i = np.arange(NB)
    di = np.abs(i - pi)
    di = np.minimum(di, NB - di)
    dj = np.abs(i - pj)
    dj = np.minimum(dj, NB - dj)
    return np.sqrt(di[:, None] ** 2 + dj[None, :] ** 2)


def main():
    os.makedirs(os.path.join(GROUND, "residual"), exist_ok=True)
    os.makedirs(os.path.join(FIN, "weights_mass"), exist_ok=True)
    os.makedirs(os.path.join(FIN, "base_mass"), exist_ok=True)

    with open(os.path.join(FIN, "dictionary.json"), encoding="utf-8") as f:
        dictionary = json.load(f)
    extras = {}
    with open(os.path.join(FIN, "contexts.jsonl"), encoding="utf-8") as f:
        for line in f:
            rec = json.loads(line)
            if rec["extra"]:
                extras[(rec["a"], rec["b"], rec["c"])] = rec["extra"]

    probs = np.load(PROBS, mmap_mode="r")
    k = ksig(NB, SIGMA)
    Kfft = np.fft.fft(k)[:, None] * np.fft.fft(k)[None, :]
    out_lines = []

    for b in range(20):
        letter = ALPHABET[b]
        peaks = dictionary[letter]
        K = len(peaks)
        res = np.zeros((20, 20, NB, NB), dtype=np.float32)
        wmass = np.zeros((20, 20, K), dtype=np.float32)
        bmass = np.zeros((20, 20), dtype=np.float32)

        for a in range(20):
            for c in range(20):
                idx = a * 400 + b * 20 + c
                emp = np.asarray(probs[idx], dtype=np.float64)
                s = emp.sum()
                if s > 0:
                    emp = emp / s

                plist = [(p["phi"], p["psi"]) for p in peaks]
                e_extra = extras.get((ALPHABET[a], letter, ALPHABET[c]), [])
                plist += [(e["phi"], e["psi"]) for e in e_extra]
                if not plist:
                    res[a, c] = emp.astype(np.float32)
                    bmass[a, c] = 1.0
                    continue

                dists = np.stack([torus_dist_grid(int((p[0] + 180) / BW) % NB,
                                                  int((p[1] + 180) / BW) % NB) for p in plist])
                nearest = np.argmin(dists, axis=0)
                mind = np.min(dists, axis=0)
                assigned = mind <= ASSIGN_R

                for kk in range(K):
                    wmass[a, c, kk] = emp[(nearest == kk) & assigned].sum()
                for kk in range(K, len(plist)):
                    m = float(emp[(nearest == kk) & assigned].sum())
                    e = e_extra[kk - K]
                    out_lines.append(json.dumps(
                        {"a": ALPHABET[a], "b": letter, "c": ALPHABET[c],
                         "phi": e["phi"], "psi": e["psi"], "weight": round(m, 6)},
                        ensure_ascii=False))

                residual = np.where(assigned, 0.0, emp)
                bmass[a, c] = float(residual.sum())
                r = np.fft.ifft2(np.fft.fft2(residual) * Kfft).real
                r = np.clip(r, 0.0, None)
                r = r / r.sum() if r.sum() > 0 else np.ones((NB, NB)) / (NB * NB)
                res[a, c] = r.astype(np.float32)

        np.save(os.path.join(GROUND, "residual", f"center_{letter}.npy"), res)
        np.save(os.path.join(FIN, "weights_mass", f"center_{letter}.npy"), wmass)
        np.save(os.path.join(FIN, "base_mass", f"center_{letter}.npy"), bmass)

    with open(os.path.join(FIN, "extras_mass.jsonl"), "w", encoding="utf-8") as f:
        f.write("\n".join(out_lines) + "\n")
    print(f"done; extras rows={len(out_lines)}")


if __name__ == "__main__":
    main()
