"""p_v2 — v2 连续二面角分布读取端（v2.1）。

接口
    p(phi, psi, a, b, c)       概率（格点质量，Σ=1）
    p_density(...)             概率密度（质量 / ΔA）
    grad(...)                  峰区内点 (p, ∂p/∂φ, ∂p/∂ψ, ∇log p)；边界带返回 None
    delta(...)                 到边界的有符号距离（度；允许 >0）
    is_forbidden(...)          格点是否禁区
    marginal(ctx, axis)        端点一维边缘（三元组中一个位置为 None）

约定（必须与链模型一致）
    1. 数值是「格点质量」：Σ_ij p = 1；密度 = p / ΔA，ΔA = (5°)² = 7.6154e-3 rad²，
       换算系数 1/ΔA = 131.31。
    2. 合法集直接用 forbidden.npy，不另设阈值 τ。
    3. 峰区梯度只在 interior（4×4 模板全在允许区）内有效；边界 2 格内用 delta。
    4. delta：允许区 > 0，禁区 < 0，零水平落在格心之间（已减半格）。

v2.1 只改读取端；final_version 的数据不变。
"""
import json
import os

import numpy as np

ROUND = r"D:\CodingPrograms\Protein Design\round"
OUT = os.path.join(ROUND, "result", "final_version")
ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
A2I = {c: i for i, c in enumerate(ALPHABET)}
NB = 72
BW = 5.0
CELL = np.deg2rad(BW) ** 2          # 7.6154e-3 rad²
INV_CELL = 1.0 / CELL               # 131.31
SQ2 = np.sqrt(2.0)
_CACHE_CAP = 400                    # 每个上下文的 interior/delta 缓存上限


def _catmull(u):
    """Catmull-Rom 基函数与其导数（u ∈ [0,1]）。"""
    u2 = u * u
    u3 = u2 * u
    w = (-0.5 * u3 + u2 - 0.5 * u,
         1.5 * u3 - 2.5 * u2 + 1.0,
         -1.5 * u3 + 2.0 * u2 + 0.5 * u,
         0.5 * u3 - 0.5 * u2)
    dw = (-1.5 * u2 + 2.0 * u - 0.5,
          4.5 * u2 - 5.0 * u,
          -4.5 * u2 + 4.0 * u + 0.5,
          1.5 * u2 - u)
    return w, dw


def _chamfer(mask):
    """到最近 True 格的环面距离（单位=格，8 邻域 chamfer，误差 ~4%）。"""
    d = np.where(mask, 0.0, np.inf)
    dirs = ((1, 0, 1.0), (-1, 0, 1.0), (0, 1, 1.0), (0, -1, 1.0),
            (1, 1, SQ2), (1, -1, SQ2), (-1, 1, SQ2), (-1, -1, SQ2))
    for _ in range(NB):
        old = d.copy()
        for di, dj, c in dirs:
            d = np.minimum(d, np.roll(np.roll(d, di, 0), dj, 1) + c)
        if np.array_equal(d, old):
            break
    return d


class V2Model:
    def __init__(self):
        self.probs = np.load(os.path.join(OUT, "probs.npy"), mmap_mode="r")
        self.forbidden = np.load(os.path.join(OUT, "forbidden.npy"), mmap_mode="r")
        with open(os.path.join(OUT, "meta.json"), encoding="utf-8") as f:
            self.meta = json.load(f)
        self._interior_cache = {}
        self._delta_cache = {}

    # ---------- 索引与坐标 ----------
    def index(self, a, b, c):
        return (A2I[a] * 20 + A2I[b]) * 20 + A2I[c]

    def _cell(self, phi, psi):
        x = (phi + 180.0) / BW
        y = (psi + 180.0) / BW
        i0 = int(np.floor(x)) % NB
        j0 = int(np.floor(y)) % NB
        return i0, j0, x - np.floor(x), y - np.floor(y)

    # ---------- 缓存 ----------
    def _interior(self, t):
        m = self._interior_cache.get(t)
        if m is None:
            allowed = ~np.asarray(self.forbidden[t], bool)
            m = allowed.copy()
            for di in (-2, -1, 0, 1):
                for dj in (-2, -1, 0, 1):
                    m &= np.roll(np.roll(allowed, di, 0), dj, 1)
            if len(self._interior_cache) > _CACHE_CAP:
                self._interior_cache.clear()
            self._interior_cache[t] = m
        return m

    def _delta_grid(self, t):
        d = self._delta_cache.get(t)
        if d is None:
            forb = np.asarray(self.forbidden[t], bool)
            d_out = _chamfer(forb)
            d_in = _chamfer(~forb)
            d_out = np.where(np.isfinite(d_out), d_out, NB)
            d_in = np.where(np.isfinite(d_in), d_in, NB)
            d = (np.where(forb, -(d_in - 0.5), (d_out - 0.5)) * BW).astype(np.float32)
            if len(self._delta_cache) > _CACHE_CAP:
                self._delta_cache.clear()
            self._delta_cache[t] = d
        return d

    # ---------- 概率 ----------
    def _p_bilinear(self, phi, psi, t):
        i0, j0, fx, fy = self._cell(phi, psi)
        i1, j1 = (i0 + 1) % NB, (j0 + 1) % NB
        g = self.probs[t]
        return float((1 - fx) * (1 - fy) * g[i0, j0] + fx * (1 - fy) * g[i1, j0]
                     + (1 - fx) * fy * g[i0, j1] + fx * fy * g[i1, j1])

    def p(self, phi, psi, a, b, c):
        """概率（格点质量）。峰区内部用双三次，边界带退回双线性。"""
        t = self.index(a, b, c)
        i0, j0, fx, fy = self._cell(phi, psi)
        if not self._interior(t)[i0, j0]:
            return self._p_bilinear(phi, psi, t)
        wu, _ = _catmull(fx)
        wv, _ = _catmull(fy)
        P = self.probs[t][np.ix_((i0 - 1 + np.arange(4)) % NB,
                                 (j0 - 1 + np.arange(4)) % NB)]
        return float(wu @ P @ wv)

    def p_density(self, phi, psi, a, b, c):
        """概率密度（每平方弧度）。"""
        return self.p(phi, psi, a, b, c) * INV_CELL

    def grad(self, phi, psi, a, b, c):
        """峰区内点：返回 (p, ∂p/∂φ, ∂p/∂ψ, ∇log p)，单位「每度」。

        边界 2 格内返回 None —— 那里请用 delta()。
        """
        t = self.index(a, b, c)
        i0, j0, fx, fy = self._cell(phi, psi)
        if not self._interior(t)[i0, j0]:
            return None
        wu, du = _catmull(fx)
        wv, dv = _catmull(fy)
        P = self.probs[t][np.ix_((i0 - 1 + np.arange(4)) % NB,
                                 (j0 - 1 + np.arange(4)) % NB)]
        p = float(wu @ P @ wv)
        gx = float(du @ P @ wv) / BW
        gy = float(wu @ P @ dv) / BW
        return p, gx, gy, np.array([gx, gy]) / p

    # ---------- 距离与合法性 ----------
    def delta(self, phi, psi, a, b, c):
        """到允许/禁区边界的有符号距离（度）。允许区 > 0，禁区 < 0。"""
        t = self.index(a, b, c)
        d = self._delta_grid(t)
        i0, j0, fx, fy = self._cell(phi, psi)
        i1, j1 = (i0 + 1) % NB, (j0 + 1) % NB
        return float((1 - fx) * (1 - fy) * d[i0, j0] + fx * (1 - fy) * d[i1, j0]
                     + (1 - fx) * fy * d[i0, j1] + fx * fy * d[i1, j1])

    def is_forbidden(self, phi, psi, a, b, c):
        t = self.index(a, b, c)
        i0, j0, _, _ = self._cell(phi, psi)
        return bool(self.forbidden[t, i0, j0])

    # ---------- 边缘 ----------
    def marginal(self, ctx, axis, pi=None):
        """端点一维边缘。

        ctx: (a,b,c)，其中恰好一个位置为 None；按 pi 对该位置加权。
        axis=1：对 φ 求和，返回 ψ 的 72 向量；axis=0：返回 φ 的 72 向量。
        返回值是「格点质量」，Σ ≈ 1。
        """
        if pi is None:
            pi = np.full(20, 1.0 / 20.0)
        k = [i for i, x in enumerate(ctx) if x is None]
        if len(k) != 1:
            raise ValueError("ctx 必须恰好有一个 None")
        k = k[0]
        acc = np.zeros(NB, dtype=np.float64)
        for lab, w in zip(ALPHABET, pi):
            trip = list(ctx)
            trip[k] = lab
            g = np.asarray(self.probs[self.index(*trip)])
            acc += w * g.sum(axis=0 if axis == 1 else 1)
        return acc


# ---------- 自检 ----------
def _selftest():
    m = V2Model()
    a = ("A", "A", "A")
    print("== 基本取值 ==")
    for phi, psi in [(-63, -42), (65, 25), (-70, 140)]:
        print(f"p({phi},{psi}|AAA) = {m.p(phi, psi, *a):.4e}  "
              f"density = {m.p_density(phi, psi, *a):.4e}  "
              f"delta = {m.delta(phi, psi, *a):+.1f}  "
              f"forbidden = {m.is_forbidden(phi, psi, *a)}")

    print("== 梯度 vs 有限差分 ==")
    phi, psi = -63.0, -42.0
    g = m.grad(phi, psi, *a)
    if g is None:
        t = m.index(*a)
        idx = np.argwhere(m._interior(t))
        i, j = idx[len(idx) // 2]
        phi, psi = -180 + BW * i, -180 + BW * j
        g = m.grad(phi, psi, *a)
    p0, gx, gy, glog = g
    h = 0.05
    fd_x = (m.p(phi + h, psi, *a) - m.p(phi - h, psi, *a)) / (2 * h)
    fd_y = (m.p(phi, psi + h, *a) - m.p(phi, psi - h, *a)) / (2 * h)
    print(f"p={p0:.4e}  gx={gx:+.4e}  fd_x={fd_x:+.4e}  rel={abs(gx-fd_x)/max(abs(gx),1e-30):.2e}")
    print(f"             gy={gy:+.4e}  fd_y={fd_y:+.4e}  rel={abs(gy-fd_y)/max(abs(gy),1e-30):.2e}")

    print("== delta 符号与范围 ==")
    t = m.index(*a)
    forb = np.asarray(m.forbidden[t], bool)
    ii, jj = np.argwhere(forb)[0]
    pf = (-180 + BW * ii, -180 + BW * jj)
    ii, jj = np.argwhere(~forb)[0]
    pa = (-180 + BW * ii, -180 + BW * jj)
    print(f"禁区点 delta={m.delta(*pf, *a):+.1f}  允许点 delta={m.delta(*pa, *a):+.1f}")
    print(f"delta 范围: [{m.delta(-180, -180, *a):+.1f}, {m.delta(0, 0, *a):+.1f}]")

    print("== 边缘 ==")
    mp = m.marginal((None, "A", "A"), axis=1)
    mq = m.marginal(("A", "A", None), axis=0)
    print(f"ΣP(ψ|·,A,A)={mp.sum():.6f}  ΣP(φ|A,A,·)={mq.sum():.6f}")


if __name__ == "__main__":
    _selftest()
