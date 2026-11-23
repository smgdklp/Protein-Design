"""build_v2 — 连续空间二面角条件分布 p(phi,psi | m_{i-1},m_i,m_{i+1})

思路（对应需求）：
  0. 简单禁区和峰区：峰区=平滑连续密度；禁区=低于阈值处，取小正地板。
  1. 循环且平滑：原始 cos/sin -> atan2 连续角，双线性撒点到环面栅格（不是计数直方图），
     再 FFT 圆周高斯卷积（周期天然、处处可导）。
  2. 抛弃离群小成分：平滑后做连通分量，质量 < 0.5% 的孤立小岛（刻意维护的特殊构型）丢弃。
  3. 任意点 >0：p = (1-eps)*峰区 + eps*均匀地板；禁区取值为 eps/(72*72)（非 0）。
  4. 数据友好：probs.npy(float32, mmap) + forbidden.npy(uint8) + peaks.jsonl + meta.json。

输出：result/final_version/
  probs.npy      (8000,72,72) float32 归一密度，idx = (a*20+b)*20+c
  forbidden.npy  (8000,72,72) uint8   1=禁区
  peaks.jsonl    每上下文特征点 (phi,psi)
  meta.json      参数 + 汇总
  summary.txt    汇总表
"""
import json
import os
import time

import numpy as np

ROUND = r"D:\CodingPrograms\Protein Design\round"
SEG = os.path.join(ROUND, "datebase", "round_bin", "segments.jsonl")
OUT = os.path.join(ROUND, "result", "final_version")
IMG = os.path.join(OUT, "imgs")

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
A2I = {c: i for i, c in enumerate(ALPHABET)}
NB = 72
BW = 5.0
NT = 8000

SIGMA = 2.5          # 圆周高斯平滑，单位=格(5°)
FLOOR_EPS = 0.01     # 分配给禁区地板的总质量
DROP_MASS = 0.005    # 连通分量质量 < 0.5% -> 丢弃（离群小成分）
DROP_CAP = 0.05      # 丢弃总质量上限（安全阀）
FORB_REL = 0.02      # 禁区阈值 = 峰高的 2%（相对阈值，不随 n 爆炸）
THR_ABS = 1e-6
PEAK_MIN_FRAC = 0.012  # 特征点：局部极大需 >= 1.2% 峰高（尽量抓到所有弱特征）
NMS_DIST = 2.0         # 非极大抑制半径（格），12° 内的极大合并为一个特征
MAX_PEAKS = 16


def log(*a):
    print("[build_v2]", *a, flush=True)


def ksig(n, s):
    i = np.arange(n, dtype=np.float64)
    d = np.minimum(i, n - i)
    k = np.exp(-0.5 * (d / s) ** 2)
    return k / k.sum()


def parse_all():
    """解析 segments.jsonl -> 全量 (t, x, y)，x,y 是 [0,72) 的连续格坐标。"""
    t_list, x_list, y_list = [], [], []
    nseg = 0
    nres = 0
    bad = 0
    t0 = time.time()
    with open(SEG, encoding="utf-8") as f:
        for ln, line in enumerate(f):
            if not line.strip():
                continue
            try:
                rec = json.loads(line)
            except Exception:
                bad += 1
                continue
            seq = rec.get("seq")
            if not seq:
                bad += 1
                continue
            n = len(seq)
            nseg += 1
            nres += n
            if n < 3:
                continue
            cp = np.asarray(rec["cos_phi"], dtype=np.float64)
            sp = np.asarray(rec["sin_phi"], dtype=np.float64)
            cq = np.asarray(rec["cos_psi"], dtype=np.float64)
            sq = np.asarray(rec["sin_psi"], dtype=np.float64)
            codes = np.fromiter((A2I.get(c, -1) for c in seq), dtype=np.int32)
            phi = np.degrees(np.arctan2(sp, cp))[1:n - 1]
            psi = np.degrees(np.arctan2(sq, cq))[1:n - 1]
            a = codes[0:n - 2]
            b = codes[1:n - 1]
            c = codes[2:n]
            if (a < 0).any() or (b < 0).any() or (c < 0).any():
                bad += 1
                continue
            t = (a * 20 + b) * 20 + c
            t_list.append(t)
            x_list.append(((phi + 180.0) / BW).astype(np.float32))
            y_list.append(((psi + 180.0) / BW).astype(np.float32))
            if ln % 10000 == 0 and ln:
                log("parse line", ln, "samples", sum(len(x) for x in x_list))
    t = np.concatenate(t_list).astype(np.int32)
    x = np.concatenate(x_list).astype(np.float32)
    y = np.concatenate(y_list).astype(np.float32)
    log(f"parse done: segments={nseg} residues={nres} bad={bad} samples={t.size} "
        f"({time.time()-t0:.1f}s)")
    return t, x, y


def splat(t, x, y):
    """双线性撒点到 72x72 环面栅格（首阶分箱，非普通计数）。"""
    i0 = np.floor(x).astype(np.int32) % NB
    j0 = np.floor(y).astype(np.int32) % NB
    fx = (x - np.floor(x)).astype(np.float64)
    fy = (y - np.floor(y)).astype(np.float64)
    i1 = (i0 + 1) % NB
    j1 = (j0 + 1) % NB
    w00 = (1 - fx) * (1 - fy)
    w10 = fx * (1 - fy)
    w01 = (1 - fx) * fy
    w11 = fx * fy
    base = t.astype(np.int64) * (NB * NB)
    sz = NT * NB * NB
    hist = np.bincount(base + i0 * NB + j0, weights=w00, minlength=sz)
    hist += np.bincount(base + i1 * NB + j0, weights=w10, minlength=sz)
    hist += np.bincount(base + i0 * NB + j1, weights=w01, minlength=sz)
    hist += np.bincount(base + i1 * NB + j1, weights=w11, minlength=sz)
    return hist.reshape(NT, NB, NB)


def torus_components(mask):
    seen = np.zeros_like(mask, dtype=bool)
    comps = []
    for i, j in zip(*np.where(mask)):
        if seen[i, j]:
            continue
        stack = [(i, j)]
        seen[i, j] = True
        cells = []
        while stack:
            x, y = stack.pop()
            cells.append((x, y))
            for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                nx, ny = (x + dx) % NB, (y + dy) % NB
                if mask[nx, ny] and not seen[nx, ny]:
                    seen[nx, ny] = True
                    stack.append((nx, ny))
        comps.append(cells)
    return comps


def local_maxima(a):
    lm = np.ones_like(a, dtype=bool)
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            if dx == 0 and dy == 0:
                continue
            lm &= a > np.roll(np.roll(a, dx, axis=0), dy, axis=1)
    return lm


def main():
    os.makedirs(IMG, exist_ok=True)
    t, x, y = parse_all()

    log("splat ...")
    t0 = time.time()
    hist = splat(t, x, y)
    del t, x, y
    log(f"splat done ({time.time()-t0:.1f}s)")

    k = ksig(NB, SIGMA)
    K = np.fft.fft(k)[:, None] * np.fft.fft(k)[None, :]

    n = hist.sum(axis=(1, 2)).astype(np.float64)

    probs = np.empty((NT, NB, NB), dtype=np.float32)
    forbidden = np.zeros((NT, NB, NB), dtype=np.uint8)
    peaks_out = []
    dropped_mass_all = []
    forb_frac_all = []
    npeaks_all = []

    CHUNK = 256
    for c0 in range(0, NT, CHUNK):
        c1 = min(c0 + CHUNK, NT)
        G = np.fft.ifft2(np.fft.fft2(hist[c0:c1].astype(np.float64), axes=(1, 2)) * K,
                         axes=(1, 2)).real
        G = np.clip(G, 0.0, None)
        for idx in range(c0, c1):
            g = G[idx - c0]
            nn = n[idx]
            s = g.sum()
            g = g / s if s > 0 else g

            # 1) 禁区阈值（=峰高的 2%，相对阈值），连通分量
            thr = max(FORB_REL * g.max(), THR_ABS)
            occ = g > thr
            p_peak = g.copy()
            dropped = 0.0
            if occ.any():
                comps = torus_components(occ)
                masses = [p_peak[tuple(zip(*comp))].sum() for comp in comps]
                order = np.argsort(masses)
                for oi in order:
                    m = masses[oi]
                    if m < DROP_MASS and dropped + m < DROP_CAP:
                        comp = comps[oi]
                        p_peak[tuple(zip(*comp))] = 0.0
                        dropped += m

            # 2) 禁区掩码 = 峰区里低于阈值的格子（含被丢弃分量）
            forb = p_peak < thr
            p_peak[forb] = 0.0
            forb_arr = forb.astype(np.uint8)

            # 3) 峰区归一 + 地板（峰区为空则退化为均匀，保证 Σ=1）
            s2 = p_peak.sum()
            floor = FLOOR_EPS / (NB * NB)
            if s2 > 0:
                p_peak = p_peak / s2
                p = (1.0 - FLOOR_EPS) * p_peak + floor
            else:
                p = np.full((NB, NB), 1.0 / (NB * NB), dtype=np.float64)

            # 4) 特征点：所有显著局部极大（1.2% 峰高）+ 非极大抑制
            peaks = []
            if p_peak.max() > 0:
                lm = local_maxima(p_peak)
                thr_peak = PEAK_MIN_FRAC * p_peak.max()
                cand = np.argwhere(lm & (p_peak > thr_peak))
                if cand.size:
                    vals = p_peak[cand[:, 0], cand[:, 1]]
                    order = np.argsort(-vals)
                    kept = []
                    for oi in order:
                        i, j = int(cand[oi, 0]), int(cand[oi, 1])
                        if p_peak[i, j] < thr_peak:
                            break
                        ok = True
                        for (ki, kj) in kept:
                            di = min(abs(i - ki), NB - abs(i - ki))
                            dj = min(abs(j - kj), NB - abs(j - kj))
                            if di * di + dj * dj <= NMS_DIST * NMS_DIST:
                                ok = False
                                break
                        if ok:
                            kept.append((i, j))
                            if len(kept) >= MAX_PEAKS:
                                break
                    peaks = [(-180.0 + BW * i + BW / 2, -180.0 + BW * j + BW / 2)
                             for i, j in kept]

            probs[idx] = p.astype(np.float32)
            forbidden[idx] = forb_arr
            peaks_out.append(peaks)
            dropped_mass_all.append(dropped)
            forb_frac_all.append(float(forb_arr.mean()))
            npeaks_all.append(len(peaks))

        log(f"modeled {c1}/{NT}")

    log("saving ...")
    np.save(os.path.join(OUT, "probs.npy"), probs)
    np.save(os.path.join(OUT, "forbidden.npy"), forbidden)

    with open(os.path.join(OUT, "peaks.jsonl"), "w", encoding="utf-8") as f:
        for idx in range(NT):
            a, b, c = idx // 400, (idx // 20) % 20, idx % 20
            rec = {"i": idx, "a": ALPHABET[a], "b": ALPHABET[b], "c": ALPHABET[c],
                   "n": int(n[idx]), "peaks": [[round(p, 1), round(q, 1)]
                                               for p, q in peaks_out[idx]]}
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")

    meta = {
        "version": "v2",
        "nb": NB, "bw_deg": BW, "sigma_bins": SIGMA,
        "floor_eps": FLOOR_EPS, "floor_value": FLOOR_EPS / (NB * NB),
        "drop_component_mass": DROP_MASS, "drop_cap": DROP_CAP,
        "forbidden_rel_threshold": FORB_REL,
        "peak_min_frac": PEAK_MIN_FRAC, "nms_dist": NMS_DIST, "max_peaks": MAX_PEAKS,
        "alphabet": ALPHABET, "ordering": "idx = (a*20+b)*20+c",
        "total_samples": int(n.sum()),
        "files": ["probs.npy", "forbidden.npy", "peaks.jsonl", "meta.json"],
    }
    with open(os.path.join(OUT, "meta.json"), "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)

    # 汇总
    npos = n[n > 0]
    dm = np.array(dropped_mass_all)
    ff = np.array(forb_frac_all)
    npk = np.array(npeaks_all)
    lines = [
        "v2 continuous dihedral distribution",
        f"sigma={SIGMA} bins floor_eps={FLOOR_EPS} drop_mass<{DROP_MASS} forb_rel={FORB_REL}",
        f"total samples = {int(n.sum())}",
        f"n per context: min={npos.min():.0f} p10={np.percentile(npos,10):.0f} "
        f"med={np.median(npos):.0f} p90={np.percentile(npos,90):.0f} max={npos.max():.0f}",
        f"dropped outlier mass per ctx: mean={dm.mean():.4f} max={dm.max():.4f}",
        f"forbidden fraction per ctx: mean={ff.mean():.3f} med={np.median(ff):.3f} "
        f"p90={np.percentile(ff,90):.3f}",
        f"feature points per ctx: mean={npk.mean():.2f} med={np.median(npk):.0f} "
        f"max={npk.max():.0f}",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT, "summary.txt"), "w", encoding="utf-8") as f:
        f.write(text + "\n")
    log("done")


if __name__ == "__main__":
    main()
