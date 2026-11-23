"""Contact sheets: all 400 context bases per center residue (20x20 mosaic)."""
import os

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROUND = r"D:\CodingPrograms\Protein Design\round"
GROUND = os.path.join(ROUND, "result", "ground", "residual_v4")
OUT = os.path.join(ROUND, "result", "fin", "base_sheets")
ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
NB = 72
GAP = 2
STEP = NB + GAP


def main():
    os.makedirs(OUT, exist_ok=True)
    fig_w = 20 * STEP / 100.0
    for ch in ALPHABET:
        res = np.load(os.path.join(GROUND, f"center_{ch}.npy"))
        big = np.zeros((20 * STEP, 20 * STEP), dtype=np.float32)
        for a in range(20):
            for c in range(20):
                t = np.asarray(res[a, c], dtype=np.float32)
                mx = t.max()
                if mx > 0:
                    t = t / mx
                big[a * STEP:a * STEP + NB, c * STEP:c * STEP + NB] = t
        fig, ax = plt.subplots(figsize=(fig_w, fig_w), dpi=100)
        ax.imshow(big, cmap="magma", interpolation="nearest")
        ax.set_title(f"center {ch}: 20x20 contexts (rows=m-1, cols=m+1), each tile normalised",
                     fontsize=10)
        ax.set_xticks(np.arange(20) * STEP + NB / 2)
        ax.set_xticklabels(list(ALPHABET), fontsize=6)
        ax.set_yticks(np.arange(20) * STEP + NB / 2)
        ax.set_yticklabels(list(ALPHABET), fontsize=6)
        fig.tight_layout()
        fig.savefig(os.path.join(OUT, f"center_{ch}.png"))
        plt.close(fig)
        print("sheet", ch)


if __name__ == "__main__":
    main()


