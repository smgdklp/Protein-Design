"""plot_v2 — 随机 500 张评估图：热图(峰区) + 禁区涂色 + 特征点标记。

输出 result/final_version/imgs/<a>-<b>-<c>.png
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
OUT = os.path.join(ROUND, "result", "final_version")
IMG = os.path.join(OUT, "imgs")
ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
NB = 72
BW = 5.0
N = 500
SEED = 1234

FORB_CMAP = ListedColormap(["#00000000", "#12121e"])


def main():
    probs = np.load(os.path.join(OUT, "probs.npy"), mmap_mode="r")
    forb = np.load(os.path.join(OUT, "forbidden.npy"), mmap_mode="r")
    with open(os.path.join(OUT, "peaks.jsonl"), encoding="utf-8") as f:
        rec_by_idx = {json.loads(l)["i"]: json.loads(l) for l in f}

    rng = np.random.default_rng(SEED)
    idxs = rng.choice(len(probs), size=N, replace=False)
    idxs = np.sort(idxs)

    extent = [-180, 180, -180, 180]
    for k, idx in enumerate(idxs):
        a, b, c = idx // 400, (idx // 20) % 20, idx % 20
        p = np.asarray(probs[idx], dtype=np.float64)
        fm = np.asarray(forb[idx], dtype=bool)
        peaks = rec_by_idx.get(int(idx), {}).get("peaks", [])
        nn = rec_by_idx.get(int(idx), {}).get("n", "")

        z = np.log10(np.clip(p, 1e-12, None))
        # probs 是 [phi, psi]；imshow 行=y、列=x，转置后 x=phi, y=psi，与标记一致
        shown = np.where(fm, np.nan, z).T

        fig, ax = plt.subplots(figsize=(3.4, 3.4), dpi=110)
        ax.imshow(shown, extent=extent, origin="lower", cmap="turbo",
                  vmin=-6, vmax=-1.3, interpolation="bilinear", aspect="equal")
        ax.imshow(fm.T.astype(float), extent=extent, origin="lower", cmap=FORB_CMAP,
                  vmin=0, vmax=1, interpolation="nearest", aspect="equal")
        for (phi, psi) in peaks:
            ax.plot(phi, psi, marker="+", color="#ff2b2b", markersize=9,
                    markeredgewidth=1.3, linestyle="none")
        ax.set_xlim(-180, 180)
        ax.set_ylim(-180, 180)
        ax.set_xlabel("phi (deg)", fontsize=7)
        ax.set_ylabel("psi (deg)", fontsize=7)
        ax.set_title(f"{ALPHABET[a]}-{ALPHABET[b]}-{ALPHABET[c]}   "
                     f"n={nn}   peaks={len(peaks)}   forbid={fm.mean()*100:.0f}%", fontsize=8)
        ax.tick_params(labelsize=6)
        fig.tight_layout()
        path = os.path.join(IMG, f"{ALPHABET[a]}-{ALPHABET[b]}-{ALPHABET[c]}.png")
        tmp = path + ".tmp"
        ok = False
        for _ in range(4):
            try:
                fig.savefig(tmp, format="png")
                os.replace(tmp, path)
                ok = True
                break
            except OSError:
                time.sleep(0.15)
        if not ok:
            print("SAVE FAILED", path, flush=True)
        plt.close(fig)
        if (k + 1) % 50 == 0:
            print(f"rendered {k+1}/{N}", flush=True)
    print("done", N, "images", flush=True)


if __name__ == "__main__":
    main()
