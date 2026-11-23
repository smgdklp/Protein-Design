"""Pooled region masses for selected neighbour contexts.

Alphabet order: A C D E F G H I K L M N P Q R S T V W Y
G=5, I=7, P=12, T=16, V=17
"""
import numpy as np

NB = 72
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
            m = np.clip(m, 1e-300, None)
            S[a, b, c] = m / m.sum()

i = np.arange(NB)[:, None]
j = np.arange(NB)[None, :]
regions = {
    "alpha": (i >= 16) & (i < 30) & (j >= 22) & (j < 38),
    "beta": (i < 18) & (j >= 54),
    "ppii": (i >= 18) & (i < 28) & (j >= 60),
    "left": (i >= 42) & (i < 58) & (j >= 32) & (j < 40),
}


def masses(p):
    return {name: float(p[m].sum()) for name, m in regions.items()}


def pooled(sel):
    num = np.zeros((NB, NB))
    den = 0.0
    for (a, b, c) in sel:
        w = occ[a, b, c]
        if w <= 0:
            continue
        num += w * S[a, b, c]
        den += w
    return num / max(den, 1e-12)


ALL = [(a, b, c) for a in range(20) for b in range(20) for c in range(20)]
print("overall", {k2: round(v, 3) for k2, v in masses(pooled(ALL)).items()})
for name, sel in [
    ("rightP", [(a, b, 12) for a in range(20) for b in range(20)]),
    ("leftP", [(12, b, c) for b in range(20) for c in range(20)]),
    ("rightG", [(a, b, 5) for a in range(20) for b in range(20)]),
    ("leftG", [(5, b, c) for b in range(20) for c in range(20)]),
    ("leftVIT", [(a, b, c) for a in (17, 7, 16) for b in range(20) for c in range(20)]),
    ("centerP", [(a, 12, c) for a in range(20) for c in range(20)]),
    ("centerG", [(a, 5, c) for a in range(20) for c in range(20)]),
]:
    p = pooled(sel)
    print(name, {k2: round(v, 3) for k2, v in masses(p).items()})
