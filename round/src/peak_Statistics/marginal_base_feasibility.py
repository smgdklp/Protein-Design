"""Feasibility of using the neighbour-agnostic marginal as the base.

For each central residue:
  - marginal support = cells covering 99% of the marginal mass
  - for every context: fraction of context mass outside the marginal support,
    and IoU between context and marginal 99% supports
Also reports how much the low-density base differs across contexts.
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


k = ksig(NB, 1.5)
K = np.fft.fft(k)[:, None] * np.fft.fft(k)[None, :]
S = np.empty_like(R)
for a in range(20):
    for b in range(20):
        for c in range(20):
            m = np.fft.ifft2(np.fft.fft2(R[a, b, c]) * K).real
            m = np.clip(m, 0.0, None)
            S[a, b, c] = m / m.sum() if m.sum() > 0 else m


def support99(m):
    v = np.sort(m.ravel())[::-1]
    c = np.cumsum(v)
    n = int(np.searchsorted(c, 0.99) + 1)
    thr = v[n - 1]
    return m >= thr


print(f"{'ctr':>3} {'outside_mass':>12} {'IoU':>6} {'marg_support_cells':>17}")
for b in range(20):
    marg = np.zeros((NB, NB))
    wsum = 0.0
    for a in range(20):
        for c in range(20):
            if occ[a, b, c] > 0:
                marg += occ[a, b, c] * S[a, b, c]
                wsum += occ[a, b, c]
    marg /= max(wsum, 1e-12)
    ms = support99(marg)
    outside = []
    ious = []
    for a in range(20):
        for c in range(20):
            if occ[a, b, c] <= 0:
                continue
            p = S[a, b, c]
            outside.append(float(p[~ms].sum()))
            ps = support99(p)
            inter = np.logical_and(ps, ms).sum()
            union = np.logical_or(ps, ms).sum()
            ious.append(inter / max(union, 1))
    print(f"{ALPHABET[b]:>3} {np.mean(outside):>12.3f} {np.mean(ious):>6.2f} {int(ms.sum()):>17}")
