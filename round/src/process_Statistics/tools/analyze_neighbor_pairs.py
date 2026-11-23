"""相邻残基对中心残基 (φ,ψ) 分布的影响分析。

数据源：result/static/triplet_probs.npy（8000×72×72 经验直方图）。
不依赖样本计数文件：用"非零格数"作为每个上下文样本量的代理权重，
并过滤掉非零格过少的上下文（噪声大）。

输出：result/static/neighbor_analysis.txt
"""

import numpy as np

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
NB = 72
BW = 5.0
SIGMA = 1.5
MIN_OCC = 100

PROBS = r"D:\CodingPrograms\Protein Design\round\result\static\triplet_probs.npy"
OUT = r"D:\CodingPrograms\Protein Design\round\result\static\neighbor_analysis.txt"


def kernel(n, sigma):
    idx = np.arange(n, dtype=float)
    d = np.minimum(idx, n - idx)
    k = np.exp(-0.5 * (d / sigma) ** 2)
    return k / k.sum()


def smooth(mat, sigma=SIGMA):
    k = kernel(NB, sigma)
    K = np.fft.fft(k)[:, None] * np.fft.fft(k)[None, :]
    m = np.fft.ifft2(np.fft.fft2(mat) * K).real
    np.clip(m, 1e-300, None, out=m)
    return m


def normalize(m):
    s = m.sum()
    return m / s if s > 0 else m


def kl(p, q):
    p = p.ravel()
    q = q.ravel()
    mask = p > 0
    return float(np.sum(p[mask] * np.log(p[mask] / q[mask])))


def js(p, q):
    m = 0.5 * (p + q)
    return 0.5 * kl(p, m) + 0.5 * kl(q, m)


def entropy(p):
    p = p.ravel()
    p = p[p > 0]
    return float(-np.sum(p * np.log(p)))


def region_masses(p):
    i = np.arange(NB)[:, None]
    j = np.arange(NB)[None, :]
    alpha = (i >= 16) & (i < 30) & (j >= 22) & (j < 38)
    beta = (i >= 0) & (i < 18) & (j >= 54) & (j < 72)
    ppii = (i >= 18) & (i < 28) & (j >= 60) & (j < 72)
    left = (i >= 42) & (i < 58) & (j >= 32) & (j < 40)
    return {
        "alpha": float(p[alpha].sum()),
        "beta": float(p[beta].sum()),
        "ppii": float(p[ppii].sum()),
        "left": float(p[left].sum()),
    }


def main():
    probs = np.load(PROBS, mmap_mode="r")
    shape = (20, 20, 20, NB, NB)
    raw = np.asarray(probs, dtype=np.float64).reshape(shape)

    occ = (raw > 0).sum(axis=(3, 4))  # 非零格数，样本量代理
    smoothed = np.empty(shape, dtype=np.float64)
    for a in range(20):
        for b in range(20):
            for c in range(20):
                smoothed[a, b, c] = normalize(smooth(raw[a, b, c]))

    w = occ.astype(np.float64)
    w[occ < MIN_OCC] = 0.0

    # 每个中心残基的边缘分布
    marg = np.empty((20, NB, NB))
    for b in range(20):
        num = np.zeros((NB, NB))
        den = 0.0
        for a in range(20):
            for c in range(20):
                ww = w[a, b, c]
                if ww <= 0:
                    continue
                num += ww * smoothed[a, b, c]
                den += ww
        marg[b] = normalize(num / max(den, 1e-12))

    lines = []
    lines.append("相邻残基对中心残基 (phi,psi) 分布的影响")
    lines.append(f"平滑 sigma={SIGMA} 格；权重=非零格数；过滤阈值 MIN_OCC={MIN_OCC}")
    lines.append("中心 灵敏度JS  左邻JS  右邻JS  边缘熵H")
    sens_all = []
    left_all = []
    right_all = []
    for b in range(20):
        # 上下文敏感度
        js_ctx = []
        ww_ctx = []
        for a in range(20):
            for c in range(20):
                ww = w[a, b, c]
                if ww <= 0:
                    continue
                js_ctx.append(js(smoothed[a, b, c], marg[b]))
                ww_ctx.append(ww)
        sens = float(np.average(js_ctx, weights=ww_ctx)) if js_ctx else 0.0

        # 左邻：固定 a，合并 c
        left_js = []
        left_w = []
        for a in range(20):
            num = np.zeros((NB, NB))
            den = 0.0
            for c in range(20):
                ww = w[a, b, c]
                if ww <= 0:
                    continue
                num += ww * smoothed[a, b, c]
                den += ww
            if den <= 0:
                continue
            left_js.append(js(normalize(num / den), marg[b]))
            left_w.append(den)
        ljs = float(np.average(left_js, weights=left_w)) if left_js else 0.0

        # 右邻：固定 c，合并 a
        right_js = []
        right_w = []
        for c in range(20):
            num = np.zeros((NB, NB))
            den = 0.0
            for a in range(20):
                ww = w[a, b, c]
                if ww <= 0:
                    continue
                num += ww * smoothed[a, b, c]
                den += ww
            if den <= 0:
                continue
            right_js.append(js(normalize(num / den), marg[b]))
            right_w.append(den)
        rjs = float(np.average(right_js, weights=right_w)) if right_js else 0.0

        sens_all.append(sens)
        left_all.append(ljs)
        right_all.append(rjs)
        lines.append(f"{ALPHABET[b]}  {sens:.4f}  {ljs:.4f}  {rjs:.4f}  {entropy(marg[b]):.3f}")

    lines.append("")
    lines.append(f"平均: 灵敏度={np.mean(sens_all):.4f} 左邻={np.mean(left_all):.4f} 右邻={np.mean(right_all):.4f}")

    # 全局：左/右邻影响力对比（按中心加权）
    lines.append("")
    lines.append("每个中心残基：左邻/右邻最大影响的残基")
    for b in range(20):
        left_effect = {}
        for a in range(20):
            num = np.zeros((NB, NB))
            den = 0.0
            for c in range(20):
                ww = w[a, b, c]
                if ww <= 0:
                    continue
                num += ww * smoothed[a, b, c]
                den += ww
            if den > 0:
                left_effect[a] = js(normalize(num / den), marg[b])
        right_effect = {}
        for c in range(20):
            num = np.zeros((NB, NB))
            den = 0.0
            for a in range(20):
                ww = w[a, b, c]
                if ww <= 0:
                    continue
                num += ww * smoothed[a, b, c]
                den += ww
            if den > 0:
                right_effect[c] = js(normalize(num / den), marg[b])
        la = max(left_effect, key=left_effect.get) if left_effect else 0
        rc = max(right_effect, key=right_effect.get) if right_effect else 0
        lines.append(
            f"{ALPHABET[b]}: 左 {ALPHABET[la]}({left_effect.get(la, 0):.4f})  "
            f"右 {ALPHABET[rc]}({right_effect.get(rc, 0):.4f})")

    # 典型上下文：P/G/β-branched 邻居
    lines.append("")
    lines.append("典型上下文区域质量 (alpha,beta,ppii,left) 与 JS")
    examples = [
        ("P", "A"), ("P", "G"), ("P", "V"), ("G", "A"), ("A", "G"),
        ("V", "A"), ("I", "A"), ("T", "A"), ("G", "G"), ("P", "P"),
    ]
    for center_letter, neigh_letter in examples:
        b = ALPHABET.index(center_letter)
        n = ALPHABET.index(neigh_letter)
        for side in ("left", "right"):
            num = np.zeros((NB, NB))
            den = 0.0
            for other in range(20):
                a, c = (n, other) if side == "left" else (other, n)
                ww = w[a, b, c]
                if ww <= 0:
                    continue
                num += ww * smoothed[a, b, c]
                den += ww
            if den <= 0:
                continue
            p = normalize(num / den)
            m = region_masses(p)
            lines.append(
                f"中心 {center_letter}, {side} {neigh_letter}: "
                f"a={m['alpha']:.2f} b={m['beta']:.2f} p={m['ppii']:.2f} L={m['left']:.2f} "
                f"JS={js(p, marg[b]):.4f}")

    # 每个中心最"异常"的相邻对（左右各取 3）
    lines.append("")
    lines.append("每个中心影响最大的 (左邻,右邻) 组合 Top3")
    for b in range(20):
        scored = []
        for a in range(20):
            for c in range(20):
                ww = w[a, b, c]
                if ww <= 0:
                    continue
                scored.append((js(smoothed[a, b, c], marg[b]), a, c, ww))
        scored.sort(reverse=True)
        top = scored[:3]
        txt = "  ".join(
            f"{ALPHABET[a]}-{ALPHABET[b]}-{ALPHABET[c]}({s:.3f})" for s, a, c, _ in top)
        lines.append(f"{ALPHABET[b]}: {txt}")

    text = "\n".join(lines)
    print(text)
    with open(OUT, "w", encoding="utf-8") as f:
        f.write(text + "\n")


if __name__ == "__main__":
    main()
