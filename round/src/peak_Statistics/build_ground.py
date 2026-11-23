"""Build the smoothed empirical ground (base) and the fin weight tables.

ground/
  base/center_X.npy        smoothed empirical marginal (99% mass kept)
  base_nopeak/center_X.npy peak-masked base used by p(phi,psi|a,b,c)
  mask/center_X.npy        peak mask (1 = removed)
  summary.txt

fin/
  dictionary.json          copy of the peak dictionary
  weights/center_X.npy     (20,20,K) matched peak weights
  contexts.jsonl           copy (matched + extra per context)
  base_nopeak/center_X.npy copy of the base used by the function example
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
PEAKDIR = os.path.join(ROUND, "result", "peak", "fin_peaks")
GROUND = os.path.join(ROUND, "result", "ground")
FIN = os.path.join(ROUND, "result", "fin")

SIGMA = 1.5
KEEP_MASS = 0.99
PEAK_MASK_R = 3.0


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


def truncate(m, keep):
    v = np.sort(m.ravel())[::-1]
    c = np.cumsum(v)
    n = int(np.searchsorted(c, keep) + 1)
    thr = v[min(n, v.size) - 1]
    out = np.where(m >= thr, m, 0.0)
    return out


def peak_mask_for(peaks):
    mask = np.zeros((NB, NB), dtype=bool)
    for p in peaks:
        i = int((p["phi"] + 180.0) / BW) % NB
        j = int((p["psi"] + 180.0) / BW) % NB
        for di in range(-int(PEAK_MASK_R), int(PEAK_MASK_R) + 1):
            for dj in range(-int(PEAK_MASK_R), int(PEAK_MASK_R) + 1):
                if di * di + dj * dj <= PEAK_MASK_R * PEAK_MASK_R:
                    mask[(i + di) % NB, (j + dj) % NB] = True
    return mask


def main():
    for d in ("base", "base_nopeak", "mask"):
        os.makedirs(os.path.join(GROUND, d), exist_ok=True)
    os.makedirs(os.path.join(FIN, "weights"), exist_ok=True)
    os.makedirs(os.path.join(FIN, "base_nopeak"), exist_ok=True)

    with open(os.path.join(PEAKDIR, "dictionary.json"), encoding="utf-8") as f:
        dictionary = json.load(f)

    probs = np.load(PROBS, mmap_mode="r")
    R = np.asarray(probs, dtype=np.float64).reshape(20, 20, 20, NB, NB)
    occ = (R > 0).sum(axis=(3, 4)).astype(float)
    k = ksig(NB, SIGMA)
    K = np.fft.fft(k)[:, None] * np.fft.fft(k)[None, :]

    lines = [f"ground: sigma={SIGMA}, keep={KEEP_MASS}, peak_mask_r={PEAK_MASK_R} bins"]
    lines.append(f"{'ctr':>3} {'mass_after_trunc':>17} {'base_mass_nopeak':>17} {'peak_cells':>10}")

    for b in range(20):
        letter = ALPHABET[b]
        marg = np.zeros((NB, NB))
        wsum = 0.0
        for a in range(20):
            for c in range(20):
                if occ[a, b, c] > 0:
                    marg += occ[a, b, c] * R[a, b, c]
                    wsum += occ[a, b, c]
        marg /= max(wsum, 1e-12)

        kept = truncate(marg, KEEP_MASS)
        mass_after = kept.sum()
        base_full = smooth(kept, K)

        mask = peak_mask_for(dictionary[letter])
        base_np = np.where(mask, 0.0, base_full)
        s = base_np.sum()
        if s > 0:
            base_np = base_np / s
        else:
            base_np = np.ones((NB, NB)) / (NB * NB)

        np.save(os.path.join(GROUND, "base", f"center_{letter}.npy"), base_full.astype(np.float32))
        np.save(os.path.join(GROUND, "base_nopeak", f"center_{letter}.npy"), base_np.astype(np.float32))
        np.save(os.path.join(GROUND, "mask", f"center_{letter}.npy"), mask)
        np.save(os.path.join(FIN, "base_nopeak", f"center_{letter}.npy"), base_np.astype(np.float32))

        # dense matched weights (20,20,K)
        ctx = {}
        with open(os.path.join(PEAKDIR, "contexts.jsonl"), encoding="utf-8") as f:
            for line in f:
                rec = json.loads(line)
                if rec["b"] == letter:
                    ctx[(rec["a"], rec["c"])] = rec
        w = np.zeros((20, 20, len(dictionary[letter])), dtype=np.float32)
        for (a, c), rec in ctx.items():
            ai, ci = ALPHABET.index(a), ALPHABET.index(c)
            for pid, val in rec["matched"].items():
                kk = int(pid[len(letter):]) - 1
                if 0 <= kk < w.shape[2]:
                    w[ai, ci, kk] = val
        np.save(os.path.join(FIN, "weights", f"center_{letter}.npy"), w)

        lines.append(f"{letter:>3} {mass_after:>17.4f} {base_np.sum():>17.4f} {int(mask.sum()):>10}")

    shutil.copy2(os.path.join(PEAKDIR, "dictionary.json"), os.path.join(FIN, "dictionary.json"))
    shutil.copy2(os.path.join(PEAKDIR, "contexts.jsonl"), os.path.join(FIN, "contexts.jsonl"))

    text = "\n".join(lines)
    print(text)
    with open(os.path.join(GROUND, "summary.txt"), "w", encoding="utf-8") as f:
        f.write(text + "\n")


if __name__ == "__main__":
    main()
