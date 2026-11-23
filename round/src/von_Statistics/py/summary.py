# -*- coding: utf-8 -*-
"""von_dpmm / 汇总:三份方阵的一致性、邻居带来的信息量、分级质量。

输出纯文本到 stdout(重定向成 summary.txt 即可)。
"""
import os
import struct

from pathlib import Path

import numpy as np

# 结果目录按脚本位置推算:cpp_round/02_sum/py/ → cpp_round/out
# (不写死绝对路径;要覆盖就用环境变量 VD_ROOT)
ROOT = os.environ.get("VD_ROOT") or str(Path(__file__).resolve().parents[2] / "out")
BINS, NCTX, CELLS = 36, 8000, 36 * 36
AA = "ACDEFGHIKLMNPQRSTVWY"
EPS = 1e-12


def read_box(path):
    with open(path, "rb") as f:
        head = f.read(64)
        magic = head[0:4]
        cnt = np.fromfile(f, dtype="<u4", count=NCTX * CELLS).reshape(NCTX, CELLS)
        prb = np.fromfile(f, dtype="<f4", count=NCTX * CELLS).reshape(NCTX, CELLS)
    return magic.decode(), cnt, prb


def read_params(path):
    with open(path, "rb") as f:
        head = f.read(72)
        raw = np.frombuffer(f.read(NCTX * 288), dtype="<u4").reshape(NCTX, 72)
    tier = raw[:, 0]
    nobs = raw[:, 1].astype(np.int64)
    keff = raw[:, 2]
    ns = raw[:, 3]
    fl = raw[:, 8:].view("<f4")   # 第 8 个 u32 起是 64 个 float(16B 头 + 16B 摘要 + 256B 成分)
    return tier, nobs, keff, ns, fl


def kl(p, q):
    p = np.maximum(p, EPS); q = np.maximum(q, EPS)
    return float((p * np.log(p / q)).sum())


def entropy(p):
    p = np.maximum(p, EPS)
    return float(-(p * np.log(p)).sum())


def main():
    mg_s, cnt, emp = read_box(os.path.join(ROOT, "staic_box", "staic_box.bin"))
    mg_v, vcnt, von = read_box(os.path.join(ROOT, "von_box", "von_box.bin"))
    tier, nobs, keff, nsamp, comp = read_params(os.path.join(ROOT, "von", "von_params.bin"))
    pred = np.fromfile(os.path.join(ROOT, "von", "von_pred_grid.bin"),
                       dtype="<f4", offset=40, count=NCTX * CELLS).reshape(NCTX, CELLS)

    print("== 文件 ==")
    print(f"staic_box magic={mg_s} total={cnt.sum()}")
    print(f"von_box   magic={mg_v} total(von_counts)={vcnt.sum()}")
    print(f"参数: K' 均值 {keff.mean():.2f}  最大 {keff.max()}  未计算 {int((keff==0).sum())}")

    n = nobs.astype(np.float64)
    N = n.sum()
    print("\n== 分级 ==")
    for t, name in [(3, "A (N>=2000)"), (2, "B (200..1999)"), (1, "C (<200)"), (0, "D (空)")]:
        m = tier == t
        if m.any():
            print(f"  {name:14s} {int(m.sum()):5d} 个   样本 {int(n[m].sum()):9d} "
                  f"({100*n[m].sum()/N:5.2f}%)   K'均值 {keff[m].mean():.2f}")

    # ---------------- 一致性
    print("\n== von_box vs staic_box (逐上下文 TV) ==")
    tv = 0.5 * np.abs(von - emp).sum(axis=1)
    tv_pred = 0.5 * np.abs(pred - emp).sum(axis=1)
    tv_vp = 0.5 * np.abs(von - pred).sum(axis=1)
    buckets = [(0, 50), (50, 200), (200, 500), (500, 1000), (1000, 2000),
               (2000, 4000), (4000, 8000)]
    print(f"{'N 区间':>14} {'上下文':>6} {'TV(von_box,emp)':>16} {'TV(pred,emp)':>14} {'TV(box,pred)':>13}")
    for lo, hi in buckets:
        m = (nobs >= lo) & (nobs < hi)
        if m.sum() == 0:
            continue
        print(f"{lo:6d}-{hi:<7d} {int(m.sum()):6d} {tv[m].mean():16.4f} "
              f"{tv_pred[m].mean():14.4f} {tv_vp[m].mean():13.4f}")
    m = nobs >= 200
    print(f"  N>=200 加权平均: TV(von_box,emp)={np.average(tv[m], weights=n[m]):.4f}  "
          f"TV(pred,emp)={np.average(tv_pred[m], weights=n[m]):.4f}  "
          f"TV(box,pred)={np.average(tv_vp[m], weights=n[m]):.4f}")

    # ---------------- 熵
    H_emp = np.array([entropy(emp[c]) for c in range(NCTX)])
    H_von = np.array([entropy(von[c]) for c in range(NCTX)])
    H_pred = np.array([entropy(pred[c]) for c in range(NCTX)])
    print("\n== 熵(自然对数;均匀=%.3f) ==" % np.log(CELLS))
    for name, H in [("经验", H_emp), ("von_box", H_von), ("式7网格", H_pred)]:
        print(f"  {name:8s} 加权均值 {np.average(H, weights=n+1e-9):.3f}")

    # ---------------- 邻居信息量
    print("\n== 邻居带来的信息量(按样本加权) ==")
    # 单点边缘 P(phi,psi|m0)
    i0 = np.arange(NCTX) // 400          # m-1
    i1 = (np.arange(NCTX) // 20) % 20    # m0
    i2 = np.arange(NCTX) % 20            # m+1
    C = cnt.astype(np.float64)

    def agg(keys, nk):
        out = np.zeros((nk, CELLS))
        np.add.at(out, keys, C)
        s = out.sum(axis=1, keepdims=True)
        return out / np.maximum(s, 1)

    P1 = agg(i1, 20)                       # (phi,psi|m0)
    P01 = agg(i0 * 20 + i1, 400)           # (m-1,m0)
    P12 = agg(i1 * 20 + i2, 400)           # (m0,m+1)
    P0 = C.sum(axis=0); P0 /= P0.sum()

    kl_c1 = np.array([kl(emp[c], P1[i1[c]]) for c in range(NCTX)])
    kl_c01 = np.array([kl(emp[c], P01[i0[c] * 20 + i1[c]]) for c in range(NCTX)])
    kl_c12 = np.array([kl(emp[c], P12[i1[c] * 20 + i2[c]]) for c in range(NCTX)])
    w = n / N
    print(f"  KL[ P(.|abc) || P(.|b) ]            = {float((w*kl_c1).sum()):.4f} nat")
    print(f"  KL[ P(.|abc) || P(.|ab) ]  (还缺 c)  = {float((w*kl_c01).sum()):.4f} nat  -> 下一个邻居的边际贡献")
    print(f"  KL[ P(.|abc) || P(.|bc) ]  (还缺 a)  = {float((w*kl_c12).sum()):.4f} nat  -> 上一个邻居的边际贡献")
    print("  (上面两个越小,说明那一侧邻居越重要)")

    # ---------------- 抽查
    print("\n== 生物学抽查(G 与 P 的 φ 分布重心) ==")
    phi_deg = -180 + (np.arange(BINS) + 0.5) * 10
    for name, key in [("G (m0=G)", AA.index("G")), ("P (m0=P)", AA.index("P"))]:
        m = (i1 == key) & (nobs >= 200)
        if m.sum() == 0:
            continue
        P = (C[m] / C[m].sum(axis=1, keepdims=True)).mean(axis=0).reshape(BINS, BINS)
        marg = P.sum(axis=1)
        print(f"  {name:10s} 平均 φ 边缘峰位 {phi_deg[np.argmax(marg)]:+.0f}°  "
              f"正 φ 质量占比 {100*P[18:].sum():.1f}%")

    print("\n== 最像/最不像的上下文(von_box vs 经验, N>=500) ==")
    m = nobs >= 500
    idx = np.argsort(tv)
    sel = [c for c in idx if m[c]]
    def nm(c):
        return AA[i0[c]] + AA[i1[c]] + AA[i2[c]]
    print("  最一致:", ", ".join(f"{nm(c)}({tv[c]:.3f})" for c in sel[:5]))
    print("  最不一致:", ", ".join(f"{nm(c)}({tv[c]:.3f})" for c in sel[-5:]))


if __name__ == "__main__":
    main()

