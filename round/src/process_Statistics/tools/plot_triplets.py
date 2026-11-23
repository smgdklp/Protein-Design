"""周期核平滑 + 3D 直方图渲染；经验统计为 0 的格用灰色底面标出。

读取 statgen 产出的 triplet_probs.npy（shape = (8000, 72, 72), float32），
对每张 72x72 概率方阵做环面（周期）高斯核卷积，再画立体直方图。
零概率格不画柱子，用灰色平面标出"没有样本"的区域，避免整张图一个颜色。

用法示例：
  python plot_triplets.py --probs ../../result/static/triplet_probs.npy \
      --out ../../result/static/img --n 500
  python plot_triplets.py --probs ../../result/static/triplet_probs.npy \
      --out ../../result/static/img_all --all
"""

import argparse
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import cm, colors  # noqa: E402

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
NB = 72
BW = 5.0  # 每格 5 度


def circular_kernel(n, sigma):
    """一维循环高斯核（格为单位）。"""
    idx = np.arange(n, dtype=np.float64)
    d = np.minimum(idx, n - idx)
    k = np.exp(-0.5 * (d / sigma) ** 2)
    return k / k.sum()


def smooth_periodic(mat, sigma):
    """二维环面高斯卷积（可分离核 + FFT）。"""
    k = circular_kernel(NB, sigma)
    K = np.fft.fft(k)[:, None] * np.fft.fft(k)[None, :]
    out = np.fft.ifft2(np.fft.fft2(mat) * K).real
    np.clip(out, 0.0, None, out=out)
    return out


def prepare(mat, sigma):
    m = np.asarray(mat, dtype=np.float64)
    s = m.sum()
    if s > 0:
        m = m / s
    m = smooth_periodic(m, sigma)
    s = m.sum()
    if s > 0:
        m = m / s
    return m


def render_one(idx, probs, out_dir, sigma, style, dpi, zero_color, gamma):
    a = ALPHABET[idx // 400]
    b = ALPHABET[(idx // 20) % 20]
    c = ALPHABET[idx % 20]

    raw = np.asarray(probs[idx], dtype=np.float64)
    total = raw.sum()
    if total > 0:
        raw = raw / total
    mask0 = raw <= 0.0

    disp = prepare(raw, sigma)
    disp[mask0] = 0.0
    vmax = float(disp.max())
    if vmax <= 0.0:
        vmax = 1.0

    centers = -180.0 + BW * (np.arange(NB) + 0.5)
    X, Y = np.meshgrid(centers, centers, indexing="ij")

    fig = plt.figure(figsize=(6.6, 5.4), dpi=dpi)
    ax = fig.add_subplot(111, projection="3d")
    ax.set_facecolor("white")
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_facecolor((1.0, 1.0, 1.0, 1.0))
        axis.pane.set_edgecolor((0.8, 0.8, 0.8, 1.0))

    if style == "bar3d":
        flat = disp.ravel()
        vals = flat / vmax
        rgba = np.empty((flat.size, 4), dtype=float)
        rgba[:] = colors.to_rgba(zero_color)
        nz = flat > 0.0
        rgba[nz] = cm.viridis(np.power(vals[nz], gamma))
        # 全部格子放进同一个 bar3d：零统计格为零高度灰块，铺成完整灰底
        ax.bar3d(
            X.ravel(), Y.ravel(), np.zeros(flat.size),
            BW, BW, flat,
            color=rgba, shade=False, edgecolor="none", linewidth=0,
        )
    else:
        ax.plot_surface(
            X, Y, disp,
            cmap="viridis", rstride=1, cstride=1, linewidth=0,
            antialiased=True, vmin=0.0, vmax=vmax,
        )

    ax.set_xlabel("phi (deg)", labelpad=6)
    ax.set_ylabel("psi (deg)", labelpad=6)
    ax.set_zlabel("P", labelpad=6)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.set_zlim(-0.02 * vmax, vmax * 1.05)
    ax.set_xticks([-180, -90, 0, 90, 180])
    ax.set_yticks([-180, -90, 0, 90, 180])
    ax.view_init(elev=35, azim=-60)
    ax.set_title(f"{a}-{b}-{c}    P(phi, psi | triplet)", fontsize=10)
    fig.tight_layout()

    path = os.path.join(out_dir, f"{a}-{b}-{c}.png")
    fig.savefig(path)
    plt.close(fig)
    return path


_W = {}


def _init(probs_path, out_dir, sigma, style, dpi, zero_color, gamma):
    _W["probs"] = np.load(probs_path, mmap_mode="r")
    _W.update(out=out_dir, sigma=sigma, style=style, dpi=dpi,
              zero_color=zero_color, gamma=gamma)


def _work(idx):
    try:
        probs = _W["probs"]
        center = ALPHABET[(idx // 20) % 20]
        out = os.path.join(_W["out"], center)
        os.makedirs(out, exist_ok=True)
        return render_one(int(idx), probs, out, _W["sigma"], _W["style"], _W["dpi"],
                          _W["zero_color"], _W["gamma"])
    except Exception as exc:  # noqa: BLE001
        return f"ERROR {idx}: {exc}"


def triplet_index(a, b, c):
    return ALPHABET.index(a.upper()) * 400 + ALPHABET.index(b.upper()) * 20 + ALPHABET.index(c.upper())


def main():
    ap = argparse.ArgumentParser(description="周期平滑 + 3D 直方图")
    ap.add_argument("--probs", required=True, help="triplet_probs.npy 路径")
    ap.add_argument("--out", required=True, help="图片输出目录")
    ap.add_argument("--all", action="store_true", help="画全部 8000 张，按中心残基分 20 个文件夹")
    ap.add_argument("--n", type=int, default=500, help="随机抽样张数（--all 时忽略）")
    ap.add_argument("--seed", type=int, default=20261108)
    ap.add_argument("--sigma", type=float, default=2.0, help="周期高斯核 sigma（格）")
    ap.add_argument("--style", choices=["bar3d", "surface"], default="bar3d")
    ap.add_argument("--dpi", type=int, default=110)
    ap.add_argument("--zero-color", default="#d0d0d0", help="经验统计为 0 的格的颜色")
    ap.add_argument("--gamma", type=float, default=0.6, help="颜色映射的 gamma，<1 提亮低值")
    ap.add_argument("--workers", type=int, default=0, help="0 = 自动")
    ap.add_argument("--single", nargs=3, metavar=("A", "B", "C"), help="只画一个三元组")
    args = ap.parse_args()

    probs = np.load(args.probs, mmap_mode="r")
    if probs.ndim != 3 or probs.shape[1:] != (NB, NB):
        raise SystemExit(f"bad shape: {probs.shape}")
    os.makedirs(args.out, exist_ok=True)

    if args.single:
        idx = triplet_index(*args.single)
        center = ALPHABET[(idx // 20) % 20]
        out = os.path.join(args.out, center)
        os.makedirs(out, exist_ok=True)
        t0 = time.time()
        path = render_one(idx, probs, out, args.sigma, args.style, args.dpi,
                          args.zero_color, args.gamma)
        print(f"{path}  ({time.time() - t0:.2f}s)")
        return

    if args.all:
        idxs = np.arange(probs.shape[0])
    else:
        rng = np.random.default_rng(args.seed)
        idxs = rng.choice(probs.shape[0], size=min(args.n, probs.shape[0]), replace=False)

    n = len(idxs)
    workers = args.workers or min(8, os.cpu_count() or 1)
    t0 = time.time()
    done = 0
    errors = 0
    with ProcessPoolExecutor(
        max_workers=workers,
        initializer=_init,
        initargs=(args.probs, args.out, args.sigma, args.style, args.dpi,
                  args.zero_color, args.gamma),
    ) as ex:
        for result in ex.map(_work, idxs, chunksize=4):
            if isinstance(result, str) and result.startswith("ERROR"):
                errors += 1
                print(result, file=sys.stderr)
            else:
                done += 1
                if done % 200 == 0:
                    print(f"{done}/{n}  {time.time() - t0:.1f}s", flush=True)

    elapsed = time.time() - t0
    print(f"done {done}/{n} errors={errors} workers={workers} style={args.style} "
          f"sigma={args.sigma} elapsed={elapsed:.1f}s")


if __name__ == "__main__":
    main()
