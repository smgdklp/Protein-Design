"""compare_v2 — 20 张经验分布 vs v2 模型交叉对照图。

左：经验直方图（result/static/triplet_probs.npy，5°格，log 色标，空格留白）
    叠加模型特征点（红 +）与模型禁区边界（青色细线）
右：v2 模型（禁区涂黑 + 同一组特征点）
"""
import json
import os
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

ROUND = r"D:\CodingPrograms\Protein Design\round"
V2 = os.path.join(ROUND, "result", "final_version")
EMP = os.path.join(ROUND, "result", "static", "triplet_probs.npy")
OUT = os.path.join(V2, "compare")
ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
NB = 72
BW = 5.0
SEED = 2024
NSAMP = 15
CURATED = [("A", "A", "A"), ("G", "G", "G"), ("P", "P", "P"), ("L", "L", "L"), ("V", "V", "V")]
FORB_CMAP = ListedColormap(["#00000000", "#141420"])


def idx(a, b, c):
    return (ALPHABET.index(a) * 20 + ALPHABET.index(b)) * 20 + ALPHABET.index(c)


def save(fig, path):
    tmp = path + ".tmp"
    for _ in range(4):
        try:
            fig.savefig(tmp, format="png")
            os.replace(tmp, path)
            return
        except OSError:
            time.sleep(0.15)
    print("SAVE FAILED", path, flush=True)


def main():
    os.makedirs(OUT, exist_ok=True)
    emp = np.load(EMP, mmap_mode="r")
    probs = np.load(os.path.join(V2, "probs.npy"), mmap_mode="r")
    forb = np.load(os.path.join(V2, "forbidden.npy"), mmap_mode="r")
    with open(os.path.join(V2, "peaks.jsonl"), encoding="utf-8") as f:
        rec = {json.loads(l)["i"]: json.loads(l) for l in f}

    n = np.array([rec[i]["n"] for i in range(8000)])
    pool = [i for i in range(8000) if n[i] >= 300]
    rng = np.random.default_rng(SEED)
    cur_idx = [idx(*t) for t in CURATED]
    pool = [i for i in pool if i not in set(cur_idx)]
    pick = list(cur_idx) + list(rng.choice(pool, size=NSAMP, replace=False))
    pick = [int(p) for p in pick]

    extent = [-180, 180, -180, 180]
    xs = -180 + BW * np.arange(NB) + BW / 2
    ys = -180 + BW * np.arange(NB) + BW / 2
    X, Y = np.meshgrid(xs, ys)

    for t in pick:
        a, b, c = ALPHABET[t // 400], ALPHABET[(t // 20) % 20], ALPHABET[t % 20]
        e = np.asarray(emp[t], dtype=np.float64)
        p = np.asarray(probs[t], dtype=np.float64)
        fm = np.asarray(forb[t], dtype=bool)
        peaks = rec[t]["peaks"]

        ze = np.log10(np.where(e > 0, e, np.nan))
        zp = np.log10(np.clip(p, 1e-12, None))
        shown = np.where(fm, np.nan, zp)

        fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.5), dpi=110)
        ax = axes[0]
        ax.imshow(ze.T, extent=extent, origin="lower", cmap="turbo",
                  vmin=-6, vmax=-1.3, interpolation="nearest", aspect="equal")
        ax.contour(X, Y, fm.T.astype(float), levels=[0.5], colors=["#22e0ff"],
                   linewidths=0.7)
        for (phi, psi) in peaks:
            ax.plot(phi, psi, marker="+", color="#ff2b2b", markersize=8,
                    markeredgewidth=1.2, linestyle="none")
        ax.set_title(f"empirical  {a}-{b}-{c}  n={n[t]}  empty={np.mean(e <= 0):.0%}",
                     fontsize=8)
        ax.set_xlabel("phi", fontsize=7)
        ax.set_ylabel("psi", fontsize=7)
        ax.tick_params(labelsize=6)

        ax = axes[1]
        ax.imshow(shown.T, extent=extent, origin="lower", cmap="turbo",
                  vmin=-6, vmax=-1.3, interpolation="bilinear", aspect="equal")
        ax.imshow(fm.T.astype(float), extent=extent, origin="lower", cmap=FORB_CMAP,
                  vmin=0, vmax=1, interpolation="nearest", aspect="equal")
        for (phi, psi) in peaks:
            ax.plot(phi, psi, marker="+", color="#ff2b2b", markersize=8,
                    markeredgewidth=1.2, linestyle="none")
        ax.set_title(f"v2 model  forbid={fm.mean():.0%}  peaks={len(peaks)}", fontsize=8)
        ax.set_xlabel("phi", fontsize=7)
        ax.set_ylabel("psi", fontsize=7)
        ax.tick_params(labelsize=6)

        fig.tight_layout()
        save(fig, os.path.join(OUT, f"{a}-{b}-{c}.png"))
        plt.close(fig)
        print("done", f"{a}-{b}-{c}", flush=True)

    print("total", len(pick))


if __name__ == "__main__":
    main()
