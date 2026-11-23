"""How large are the peak position offsets relative to the shared dictionary?

Reuses the same peak extraction / clustering as peak_pipeline.py and
reports, per central residue, the distribution of (delta_phi, delta_psi)
for peaks matched to a shared cluster within 3 bins.
"""
import numpy as np

NB = 72
BW = 5.0
ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
PROBS = r"D:\CodingPrograms\Protein Design\round\result\static\triplet_probs.npy"
SIGMA = 1.5
PEAK_FRAC = 0.08
MAX_PEAKS = 6
CLUSTER_R = 3.0
MIN_SUPPORT = 10
MIN_MASS_FRAC = 0.02
MIN_OCC = 100


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
    di = abs(i1 - i2)
    dj = abs(j1 - j2)
    di = min(di, NB - di)
    dj = min(dj, NB - dj)
    return float(np.hypot(di, dj))


def extract(m):
    lm = local_maxima(m)
    thr = PEAK_FRAC * m.max()
    cand = np.argwhere(lm & (m >= thr))
    if len(cand) == 0:
        return []
    vals = m[cand[:, 0], cand[:, 1]]
    order = np.argsort(-vals)[:MAX_PEAKS]
    return [{"i": int(cand[k, 0]), "j": int(cand[k, 1]),
             "mass": float(m[cand[k, 0], cand[k, 1]])} for k in order]


P = np.load(PROBS, mmap_mode="r")
R = np.asarray(P, dtype=np.float64).reshape(20, 20, 20, NB, NB)
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

print(f"{'ctr':>3} {'matched':>7} {'med|dphi|':>9} {'med|dpsi|':>9} {'p90|dphi|':>9} "
      f"{'p90|dpsi|':>9} {'<=1bin':>6} {'<=2bin':>6} {'<=3bin':>6}")
for b in range(20):
    allp = []
    for a in range(20):
        for c in range(20):
            if occ[a, b, c] < MIN_OCC:
                continue
            for p in extract(S[a, b, c]):
                p.update({"a": a, "c": c})
                allp.append(p)
    if not allp:
        continue
    allp.sort(key=lambda p: -p["mass"])
    clusters = []
    for p in allp:
        best, bd = None, CLUSTER_R + 1e-9
        for cl in clusters:
            d = torus_dist(p["i"], p["j"], cl["i"], cl["j"])
            if d < bd:
                bd, best = d, cl
        if best is None:
            clusters.append({"i": p["i"], "j": p["j"], "mass": p["mass"], "support": 1})
        else:
            best["mass"] += p["mass"]
            best["support"] += 1
    total = sum(c["mass"] for c in clusters)
    clusters.sort(key=lambda c: -c["mass"])
    kept, cum = [], 0.0
    for cl in clusters:
        frac = cl["mass"] / total
        if cl["support"] < MIN_SUPPORT or frac < MIN_MASS_FRAC:
            continue
        kept.append(cl)
        cum += frac
        if cum >= 0.90:
            break

    dphi, dpsi = [], []
    for p in allp:
        best, bd = None, CLUSTER_R + 1e-9
        for cl in kept:
            d = torus_dist(p["i"], p["j"], cl["i"], cl["j"])
            if d < bd:
                bd, best = d, cl
        if best is None:
            continue
        di = abs(p["i"] - best["i"])
        di = min(di, NB - di)
        dj = abs(p["j"] - best["j"])
        dj = min(dj, NB - dj)
        dphi.append(di)
        dpsi.append(dj)
    if not dphi:
        continue
    dphi = np.array(dphi) * BW
    dpsi = np.array(dpsi) * BW
    print(f"{ALPHABET[b]:>3} {len(dphi):>7} {np.median(dphi):>9.1f} {np.median(dpsi):>9.1f} "
          f"{np.percentile(dphi,90):>9.1f} {np.percentile(dpsi,90):>9.1f} "
          f"{(dphi <= 5).mean():>6.2f} {(dphi <= 10).mean():>6.2f} {(dphi <= 15).mean():>6.2f}")
