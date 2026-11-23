# -*- coding: utf-8 -*-
"""von_dpmm / 画图:三组结果各随机 5 个上下文,离散 36x36 -> 环面插值连续图。

用法:
    python plot.py [--root <cpp_round>/out] [--n 5] [--seed 20240601]
不写文件到 C 盘;输出到 <root>/img。默认 root 按脚本位置推算(02_sum/py → out)。
"""
import argparse
import os
import struct
import sys

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import zoom

BINS = 36
NCTX = 8000
CELLS = BINS * BINS
AA = "ACDEFGHIKLMNPQRSTVWY"


def ctx_name(c):
    return AA[c // 400] + AA[(c // 20) % 20] + AA[c % 20]


def read_box(path):
    with open(path, "rb") as f:
        head = f.read(64)
        magic, ver, bins, nctx, ncell, enc, rsv = struct.unpack_from("<4s6I", head, 0)
        aa = head[28:48].decode("ascii", "replace")
        total, = struct.unpack_from("<Q", head, 48)
        flags, hb = struct.unpack_from("<2I", head, 56)
        assert magic in (b"SBX1", b"VBX1"), magic
        assert bins == BINS and nctx == NCTX and ncell == CELLS
        cnt = np.fromfile(f, dtype="<u4", count=NCTX * CELLS).reshape(NCTX, BINS, BINS)
        prb = np.fromfile(f, dtype="<f4", count=NCTX * CELLS).reshape(NCTX, BINS, BINS)
    return dict(magic=magic.decode(), aa=aa, total=total, counts=cnt, probs=prb)


def read_params(path):
    with open(path, "rb") as f:
        head = f.read(72)
        magic, ver, bins, nctx, kmax, naa = struct.unpack_from("<4s5I", head, 0)
        aa = head[24:44].decode("ascii", "replace")
        total, = struct.unpack_from("<Q", head, 48)
        flags, recb = struct.unpack_from("<2I", head, 56)
        assert magic == b"VPM1" and recb == 288
        raw = np.fromfile(f, dtype="<u4", count=NCTX * recb // 4).reshape(NCTX, recb // 4)
    tier = raw[:, 0]
    nobs = raw[:, 1]
    keff = raw[:, 2]
    ns = raw[:, 3]
    fl = raw[:, 4:].view("<f4").reshape(NCTX, 68)
    comp = fl[:, 4:].reshape(NCTX, 8, 8)
    return dict(tier=tier, n=nobs, k=keff, nsamp=ns,
                kpm=fl[:, 0], kpsd=fl[:, 1], ll=fl[:, 2], bic=fl[:, 3], comp=comp)


def smooth(a, f=8):
    """环面插值:grid-wrap 的三次样条,避免 ±180 处出现假断崖。"""
    return zoom(a, f, order=3, mode="grid-wrap")


def vlim(z):
    hi = float(z.max())
    if hi <= 0:
        return 1e-9
    pos = z[z > 0]
    if pos.size > 8:
        hi = min(hi, float(np.percentile(pos, 99.5)) * 1.05)
    return max(hi, 1e-9)


def panel(ax, z, title, discrete=False):
    x = np.linspace(-180, 180, BINS + 1)
    y = np.linspace(-180, 180, BINS + 1)
    vmax = vlim(z)
    if discrete:
        ax.pcolormesh(x, y, z.T, cmap="magma", shading="flat", vmin=0, vmax=vmax)
    else:
        ax.imshow(z, origin="lower", extent=[-180, 180, -180, 180],
                  cmap="magma", interpolation="bilinear", vmin=0, vmax=vmax)
        # 只画相对峰值的等值线,并把远低于峰值的噪声场掩掉
        hi = float(z.max())
        zm = np.where(z > 0.18 * hi, z, np.nan)
        lv = [hi * f for f in (0.25, 0.45, 0.65, 0.85)]
        ax.contour(np.linspace(-180, 180, z.shape[0]),
                   np.linspace(-180, 180, z.shape[1]), zm.T,
                   levels=lv, colors="w", linewidths=0.55, alpha=0.65)
    ax.set_xticks(np.arange(-180, 181, 90))
    ax.set_yticks(np.arange(-180, 181, 90))
    ax.set_xlabel(r"$\varphi$ (deg)")
    ax.set_ylabel(r"$\psi$ (deg)")
    ax.set_title(title, fontsize=8)
    ax.set_aspect("equal")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root",
                    default=os.environ.get("VD_ROOT")
                    or str(Path(__file__).resolve().parents[2] / "out"))
    ap.add_argument("--n", type=int, default=5)
    ap.add_argument("--seed", type=int, default=20240601)
    args = ap.parse_args()

    img = os.path.join(args.root, "img")
    os.makedirs(img, exist_ok=True)

    P = read_params(os.path.join(args.root, "von", "von_params.bin"))
    S = read_box(os.path.join(args.root, "staic_box", "staic_box.bin"))
    V = read_box(os.path.join(args.root, "von_box", "von_box.bin"))
    print("staic total =", S["total"], " von_box total =", V["total"])

    rng = np.random.default_rng(args.seed)
    picks = rng.choice(NCTX, size=args.n, replace=False)

    sets = [
        ("von", "式(7) 后验预测密度", P),
        ("von_box", "von 参数模型离散采样", V),
        ("staic_box", "经验直方图(未平滑)", S),
    ]

    for tag, label, data in sets:
        for i, c in enumerate(picks):
            if tag == "von":
                grid = np.stack([P["comp"][c, :, 0]]).ravel()
                # von 用预测网格文件
                z = np.fromfile(os.path.join(args.root, "von", "von_pred_grid.bin"),
                                dtype="<f4", offset=40, count=NCTX * CELLS)
                z = z.reshape(NCTX, BINS, BINS)[c]
                kmax = int(P["k"][c])
            else:
                z = data["probs"][c]
                kmax = int(P["k"][c])
            zs = smooth(z, 8)
            fig = plt.figure(figsize=(12.6, 4.3))
            ax1 = fig.add_subplot(1, 3, 1)
            ax2 = fig.add_subplot(1, 3, 2)
            ax3 = fig.add_subplot(1, 3, 3, projection="3d")
            panel(ax1, z, f"discrete 36x36  {ctx_name(c)}  N={P['n'][c]}", discrete=True)
            panel(ax2, zs, "torus interpolation (spline, grid-wrap)")
            ng = zs.shape[0]
            X, Y = np.meshgrid(np.linspace(-180, 180, ng), np.linspace(-180, 180, ng))
            ax3.plot_surface(X, Y, zs.T, cmap="magma", linewidth=0, antialiased=True, rcount=72, ccount=72)
            ax3.set_xlabel(r"$\varphi$"); ax3.set_ylabel(r"$\psi$")
            ax3.set_title("3D", fontsize=8)
            ax3.view_init(28, -125)
            fig.suptitle(f"{tag}  |  P(phi,psi | {ctx_name(c)[0]}, {ctx_name(c)[1]}, "
                         f"{ctx_name(c)[2]})  N={int(P['n'][c])}  K'={kmax}  "
                         f"tier={int(P['tier'][c])}  max={z.max():.3f}",
                         fontsize=9)
            fig.tight_layout(rect=[0, 0, 1, 0.94])
            out = os.path.join(img, f"{tag}_{i+1}_{ctx_name(c)}.png")
            fig.savefig(out, dpi=130)
            plt.close(fig)
            print("wrote", out)

    # 对照图:同一批上下文 von_box vs staic_box + 差值
    for i, c in enumerate(picks):
        zv = V["probs"][c]
        ze = S["probs"][c]
        tv = 0.5 * np.abs(zv - ze).sum()
        fig, axs = plt.subplots(1, 3, figsize=(12.6, 4.2))
        panel(axs[0], zv, f"von_box  {ctx_name(c)}  N={P['n'][c]} K'={int(P['k'][c])}")
        panel(axs[1], ze, "staic_box (empirical)")
        d = zv - ze
        m = max(float(np.abs(d).max()), 1e-9)
        im = axs[2].imshow(d.T, origin="lower", extent=[-180, 180, -180, 180],
                           cmap="RdBu_r", vmin=-m, vmax=m)
        axs[2].set_title(f"difference  TV={tv:.3f}", fontsize=8)
        axs[2].set_xlabel(r"$\varphi$ (deg)"); axs[2].set_aspect("equal")
        fig.colorbar(im, ax=axs[2], fraction=0.046)
        fig.suptitle(f"compare  {tag if False else ''}{ctx_name(c)}", fontsize=9)
        fig.tight_layout(rect=[0, 0, 1, 0.94])
        out = os.path.join(img, f"cmp_{i+1}_{ctx_name(c)}.png")
        fig.savefig(out, dpi=130)
        plt.close(fig)
        print("wrote", out)


if __name__ == "__main__":
    sys.exit(main())



