"""Same central residue, different neighbour contexts: peak height collapse.

For each central residue b:
  - weighted marginal over 400 contexts
  - mean peak probability of the individual contexts
  - peak probability of the marginal
  - number of modes and plateau size of the marginal
  - circular spread of the per-context peak positions
"""
import numpy as np

NB = 72
ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
P = np.load(r"D:\CodingPrograms\Protein Design\round\result\static\triplet_probs.npy", mmap_mode="r")
R = np.asarray(P, dtype=np.float64).reshape(20, 20, 20, NB, NB)
occ = (R > 0).sum(axis=(3, 4)).astype(float)


def ksig(n, s):
    i = np.arange(n, dtype=float)
    d = np.minimum(i, n - i)
    k = np.exp(-0.5 * (d / s) ** 2)
    return k / k.sum()


SIGMA = 1.5
k = ksig(NB, SIGMA)
K = np.fft.fft(k)[:, None] * np.fft.fft(k)[None, :]
S = np.empty_like(R)
for a in range(20):
    for b in range(20):
        for c in range(20):
            m = np.fft.ifft2(np.fft.fft2(R[a, b, c]) * K).real
            m = np.clip(m, 0.0, None)
            S[a, b, c] = m / m.sum() if m.sum() > 0 else m


def local_max_count(M):
    lm = np.ones_like(M, dtype=bool)
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            if dx == 0 and dy == 0:
                continue
            lm &= M > np.roll(np.roll(M, dx, axis=0), dy, axis=1)
    peak = M.max()
    return int((lm & (M >= 0.05 * peak)).sum()), float(peak), int((M >= 0.9 * peak).sum())


def circ_std(angles_rad, weights):
    z = np.sum(weights * np.exp(1j * angles_rad)) / np.sum(weights)
    r = np.abs(z)
    if r <= 1e-12:
        return 180.0
    return float(np.degrees(np.sqrt(-2.0 * np.log(r))))


print(f"sigma={SIGMA}; weighted by occupied cells")
print(f"{'ctr':>4} {'ctxPeak':>8} {'margPeak':>9} {'ratio':>6} {'modes':>6} {'plateau':>8} "
      f"{'sdPhi':>7} {'sdPsi':>7}")
for b in range(20):
    w = occ[:, b, :]
    wsum = w.sum()
    marg = np.zeros((NB, NB))
    peaks = []
    phis = []
    psis = []
    ww = []
    for a in range(20):
        for c in range(20):
            wwv = w[a, c]
            if wwv <= 0:
                continue
            m = S[a, b, c]
            marg += wwv * m
            idx = np.unravel_index(np.argmax(m), m.shape)
            peaks.append(float(m[idx]))
            phis.append(np.radians(-180.0 + 5.0 * idx[0] + 2.5))
            psis.append(np.radians(-180.0 + 5.0 * idx[1] + 2.5))
            ww.append(wwv)
    marg /= wsum
    ctx_peak = float(np.average(peaks, weights=ww))
    modes, marg_peak, plateau = local_max_count(marg)
    sd_phi = circ_std(np.array(phis), np.array(ww))
    sd_psi = circ_std(np.array(psis), np.array(ww))
    print(f"{ALPHABET[b]:>4} {ctx_peak:8.4f} {marg_peak:9.4f} {ctx_peak / marg_peak:6.2f} "
          f"{modes:6d} {plateau:8d} {sd_phi:7.1f} {sd_psi:7.1f}")
