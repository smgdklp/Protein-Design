# src/peak_Statistics —— 仅峰成分的统计分析

对 `P(φ,ψ | m_{i-1}, m_i, m_{i+1})` 只保留峰成分，不拟合基底：

1. 对 8000 个上下文（σ=1.5 平滑）提峰；
2. 以每个中心残基为单位，把 400 个上下文的峰位做环面聚类，得到**共享峰字典**；
3. 每个上下文只保留“往各共享峰的权重”（400×K），不保存峰形状；
4. 用中心残基的边缘分布交叉对照，统计被叠加隐藏（swallow）的峰质量；
5. 标出无法被共享峰解释的异类上下文，并从 `result/static/img_all` 拷贝对应图像。

## 运行

```powershell
# 用 process_Statistics 的 uv venv
src\process_Statistics\.venv\Scripts\python.exe src\peak_Statistics\peak_pipeline.py
```

或经 CMake 目标：

```powershell
cmake --build build --target peak_analysis
```

## 输出（`result\peak_ground\peak`）

| 文件 | 内容 |
|---|---|
| `peaks\center_X.json` | 该中心残基的共享峰：位置、质量占比、支持上下文数、峰位散布 |
| `weights\center_X.npy` | `(20, 20, K)`，索引 `[a, c, k]`，每上下文对各共享峰的峰质量权重 |
| `marginal\center_X.npy` | 中心残基边缘分布（72×72），用于交叉对照 |
| `outliers.csv` | 共享峰解释不了的异类上下文（a,b,c,未匹配比例,样本量） |
| `images\<center>\` | 异类图像，直接从 `img_all` 拷贝，不重画（当前目录不存在） |
| `report.txt` | 汇总表 |

## 参数

| 参数 | 默认 | 含义 |
|---|---|---|
| `SIGMA` | 1.5 格 | 周期平滑宽度 |
| `PEAK_FRAC` | 0.08 | 峰高需达到该上下文峰值的比例 |
| `MAX_PEAKS` | 6 | 每上下文最多取的峰数 |
| `CLUSTER_R` | 3 格（15°） | 环面峰位聚类半径 |
| `MIN_SUPPORT` | 10 | 共享峰至少被多少上下文支持 |
| `MIN_MASS_FRAC` | 0.02 | 共享峰至少占总峰质量的比例 |
| `OUTLIER_UNMATCHED` | 0.40 | 未匹配峰质量超过该比例即判为异类 |

## 判定口径

- `context_peak_coverage`：上下文峰质量中被共享峰覆盖的比例（越高越适合“一组峰描述所有组合”）；
- `marginal_peak_mass_swallowed`：共享峰里在边缘分布中找不到对应峰的质量占比（衡量叠加后隐藏了多少峰）。

## 最终表示（`result\peak_ground\peak\fin_peaks`）

```powershell
src\process_Statistics\.venv\Scripts\python.exe src\peak_Statistics\export_fin_peaks.py
```

| 文件 | 内容 |
|---|---|
| `dictionary.json` | 每个中心残基的共享峰：`id`、固定位置 `(phi,psi)`、质量占比、支持数 |
| `contexts.jsonl` | 8000 条上下文：`matched{峰id: 权重}`、`extra[{phi,psi,weight}]`、`coverage`、`occ` |
| `outliers.jsonl` | 带显著 `extra` 的上下文（`coverage<0.75` 或单个 extra 权重 ≥0.2），约 1800 条 |

- 共享峰按 `id` 引用，位置固定在字典里；**未出现的 id 权重即为 0**，不需要单独维护；
- `extra` 是字典覆盖不到的峰，自带位置和权重——所谓“异类”就是这些 extra，不再是独立概念；
- 上下文记录自洽：`coverage = Σ matched 权重`，剩余质量即 extra。
