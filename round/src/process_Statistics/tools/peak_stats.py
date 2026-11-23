"""Peak / plateau statistics of the empirical (phi,psi) histograms.

For each of the 8000 matrices and several smoothing widths:
  - number of local maxima above 5% of the peak
  - plateau size: cells >= 90% of the peak
  - peak probability and mass in the top 1% cells
Aggregates medians and per-center-residue values.
"""
import numpy as np

NB = 72
ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
P = np.load(r"D:\CodingPrograms\Protein Design\round\result\static\triplet_probs.npy", mmap_mode="r")
R = np.asarray(P, dtype=np.float32).reshape(8000, NB, NB)


def ksig(n, s):
    i = np.arange(n, dtype=float)
    d = np.minimum(i, n - i)
    k = np.exp(-0.5 * (d / s) ** 2)
    return k / k.sum()


def smooth_sigma(sigma):
    k = ksig(NB, sigma)
    K = np.fft.fft(k)[:, None] * np.fft.fft(k)[None, :]
    out = np.empty_like(R)
    for t in range(8000):
        m = np.fft.ifft2(np.fft.fft2(R[t]) * K).real
        m = np.clip(m, 0.0, None)
        s = m.sum()
        out[t] = m / s if s > 0 else m
    return out


def local_max_count(S):
    M = np.ones_like(S, dtype=bool)
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            if dx == 0 and dy == 0:
                continue
            M &= S > np.roll(np.roll(S, dx, axis=1), dy, axis=2)
    return M


for sigma in (1.0, 1.5, 2.0, 3.0):
    S = smooth_sigma(sigma)
    peak = S.max(axis=(1, 2))
    lm = local_max_count(S)
    modes = (lm & (S >= 0.05 * peak[:, None, None])).sum(axis=(1, 2))
    plateau = (S >= 0.9 * peak[:, None, None]).sum(axis=(1, 2))
    top = np.sort(S.reshape(8000, -1), axis=1)[:, -int(0.01 * NB * NB):].sum(axis=1)

    print(f"sigma={sigma}")
    print(f"  peak prob      median={np.median(peak):.4f} p90={np.percentile(peak,90):.4f}")
    print(f"  modes(>5%)     median={np.median(modes):.1f} mean={modes.mean():.1f} p90={np.percentile(modes,90):.1f}")
    print(f"  plateau(>=90%) median={np.median(plateau):.0f} mean={plateau.mean():.1f}")
    print(f"  top1% cell mass median={np.median(top):.3f}")
    # per center residue
    m_by_center = modes.reshape(20, 20, 20).mean(axis=(1, 2))
    print("  modes by center: " + " ".join(f"{ALPHABET[b]}={m_by_center[b]:.1f}" for b in range(20)))
