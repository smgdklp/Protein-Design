"""Entropy decomposition: how much do the neighbours add beyond the central residue?

H0                    : entropy of the pooled (phi,psi) distribution
H(center)             : mean entropy given the central residue
H(center,left)        : mean entropy given (left,center)
H(center,right)       : mean entropy given (center,right)
H(center,left,right)  : mean entropy given the full triplet

I = H(center) - H(...) is the information contributed by that neighbour,
given the central residue.
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


def entropy(p):
    p = p.ravel()
    p = p[p > 0]
    return float(-np.sum(p * np.log(p)))


def smooth_all(sigma):
    k = ksig(NB, sigma)
    K = np.fft.fft(k)[:, None] * np.fft.fft(k)[None, :]
    S = np.empty_like(R)
    for a in range(20):
        for b in range(20):
            for c in range(20):
                m = np.fft.ifft2(np.fft.fft2(R[a, b, c]) * K).real
                m = np.clip(m, 0.0, None)
                S[a, b, c] = m / m.sum() if m.sum() > 0 else m
    return S


def weighted_entropy(blocks, weights):
    tot = weights.sum()
    num = 0.0
    for blk, w in zip(blocks, weights):
        if w <= 0:
            continue
        num += w * entropy(blk)
    return num / tot


for sigma in (1.0, 1.5):
    S = smooth_all(sigma)
    w = occ

    pooled = np.zeros((NB, NB))
    for a in range(20):
        for b in range(20):
            for c in range(20):
                if w[a, b, c] > 0:
                    pooled += w[a, b, c] * S[a, b, c]
    H0 = entropy(pooled / pooled.sum())

    # H(center)
    blocks, weights = [], []
    for b in range(20):
        m = np.zeros((NB, NB))
        wb = 0.0
        for a in range(20):
            for c in range(20):
                if w[a, b, c] > 0:
                    m += w[a, b, c] * S[a, b, c]
                    wb += w[a, b, c]
        if wb > 0:
            blocks.append(m / wb)
            weights.append(wb)
    Hc = weighted_entropy(blocks, np.array(weights))

    # H(center,left)
    blocks, weights = [], []
    for a in range(20):
        for b in range(20):
            m = np.zeros((NB, NB))
            wb = 0.0
            for c in range(20):
                if w[a, b, c] > 0:
                    m += w[a, b, c] * S[a, b, c]
                    wb += w[a, b, c]
            if wb > 0:
                blocks.append(m / wb)
                weights.append(wb)
    Hcl = weighted_entropy(blocks, np.array(weights))

    # H(center,right)
    blocks, weights = [], []
    for b in range(20):
        for c in range(20):
            m = np.zeros((NB, NB))
            wb = 0.0
            for a in range(20):
                if w[a, b, c] > 0:
                    m += w[a, b, c] * S[a, b, c]
                    wb += w[a, b, c]
            if wb > 0:
                blocks.append(m / wb)
                weights.append(wb)
    Hcr = weighted_entropy(blocks, np.array(weights))

    # H(triplet)
    blocks, weights = [], []
    for a in range(20):
        for b in range(20):
            for c in range(20):
                if w[a, b, c] > 0:
                    blocks.append(S[a, b, c])
                    weights.append(w[a, b, c])
    Ht = weighted_entropy(blocks, np.array(weights))

    print(f"sigma={sigma}: H0={H0:.4f} H(center)={Hc:.4f} H(center,left)={Hcl:.4f} "
          f"H(center,right)={Hcr:.4f} H(triplet)={Ht:.4f}")
    print(f"  I(center)={H0 - Hc:.4f} nats")
    print(f"  I(left|center)={Hc - Hcl:.4f}  I(right|center)={Hc - Hcr:.4f}  "
          f"I(both|center)={Hc - Ht:.4f} nats")
    print(f"  bits: left={np.log2(np.e) * (Hc - Hcl):.3f} right={np.log2(np.e) * (Hc - Hcr):.3f} "
          f"both={np.log2(np.e) * (Hc - Ht):.3f}")
    print(f"  fraction of remaining entropy: left={(Hc - Hcl) / Hc:.3%} "
          f"right={(Hc - Hcr) / Hc:.3%} both={(Hc - Ht) / Hc:.3%}")
