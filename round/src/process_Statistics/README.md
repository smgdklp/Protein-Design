# process_Statistics —— 二面角条件概率统计

以相邻三残基 `(m_{i-1}, m_i, m_{i+1})` 为条件，统计 (φ, ψ) 的 72×72 离散联合分布：

```
C[a][b][c][i][j] += 1        # a=m(i-1), b=m(i), c=m(i+1)，i=φ 箱，j=ψ 箱
P(t)[i][j] = C(t)[i][j] / Σ C(t)
```

产出 `20×20×20 = 8000` 张 `72×72` 的概率方阵（纯频数统计后逐张归一）。

## 构建

```powershell
cmake -S "D:\CodingPrograms\Protein Design\round" -B "D:\CodingPrograms\Protein Design\round\build" -G "Visual Studio 18 2026" -A x64
cmake --build "D:\CodingPrograms\Protein Design\round\build" --config Release --parallel
```

产物：`round\bin\Release\statgen.exe`，`config.json` 会自动复制到 exe 同目录。

## 配置（`config.json`，与 exe 同目录）

```json
{
    "BIN_PATH": "D:\\CodingPrograms\\Protein Design\\round\\datebase\\round_bin",
    "RESULT_PATH": "D:\\CodingPrograms\\Protein Design\\round\\result\\static"
}
```

| 键 | 含义 |
|---|---|
| `BIN_PATH` | 输入目录，读取其中的 `segments.jsonl` |
| `RESULT_PATH` | 输出目录，不存在会自动创建 |

## 输出

| 文件 | 说明 |
|---|---|
| `triplet_probs.npy` | float32，shape `(8000, 72, 72)`，每张非空方阵 Σ=1 |
| `summary.txt` | 统计摘要与自检结果 |
| `statgen.log` | 运行日志 |
| 可选 `triplet_counts.npy` | 加 `--counts` 时输出 uint32 原始计数 |

**索引约定**（固定，不额外产出索引文件）：

- 残基字母表顺序：`ACDEFGHIKLMNPQRSTVWY`（20 个单字母，互不重复）
- `idx = a*400 + b*20 + c`，其中 `a=m_{i-1}`、`b=m_i`、`c=m_{i+1}`
- 矩阵第 0 维是 φ、第 1 维是 ψ

**分箱约定**：

- 72 格 × 5°，范围 `[-180°, 180°)`，第 k 格为 `[-180+5k, -180+5(k+1))`
- `bin = floor((φ+180)/5)`；恰好 +180° 回绕到 0 号箱，−180° 也落在 0 号箱
- 落在格边界上的角度一律归入右侧格，不会丢点

**统计口径**：只统计同时有 φ 和 ψ 的位点，即每段 `i = 1 … n-2`，
总样本数 = 残基总数 − 2×片段数。上游的 MSE→MET、非标准截断、最小长度 3 全部继承。

## 运行

```powershell
.\bin\Release\statgen.exe                 # 读 exe 同目录 config.json
.\bin\Release\statgen.exe --counts        # 额外写计数张量
.\bin\Release\statgen.exe --inspect A A A # 打印某个三元组的 72x72 归一矩阵(CSV)
.\bin\Release\statgen.exe --help
```

## 绘图工具（uv 托管）

```powershell
# 1) 检测 uv 工具链（首次加 -Install 会从 uv 缓存离线装 numpy/matplotlib）
powershell -ExecutionPolicy Bypass -File src\process_Statistics\tools\check_env.ps1 -Install

# 2) 周期核平滑 + 3D 直方图；随机抽 500 张
src\process_Statistics\.venv\Scripts\python.exe src\process_Statistics\tools\plot_triplets.py `
    --probs result\static\triplet_probs.npy --out result\static\img --n 500 --style bar3d --sigma 2
```

- `--sigma`：环面高斯核宽度（单位：格，默认 2 格 = 10°），平滑只用于绘图，不改变产出的概率方阵
- `--style`：`bar3d`（立体直方图，默认）或 `surface`（3D 曲面）
- `--single A A A`：只画一个三元组，用于调试
- 抽样固定随机种子（默认 20261108），图像命名 `<m-1>-<m>-<m+1>.png`

### 全量 8000 张（按中心残基分 20 个文件夹）

```powershell
src\process_Statistics\.venv\Scripts\python.exe src\process_Statistics\tools\plot_triplets.py `
    --probs result\static\triplet_probs.npy --out result\static\img_all --all
```

- 输出 `img_all\<中心残基>\<a>-<b>-<c>.png`，中心残基用单字母（A/C/D/.../Y），每个文件夹 400 张
- `--zero-color`（默认 `#d0d0d0`）：**经验统计为 0 的格**画成零高度灰块，铺成完整灰底；
  有数据的格是彩色柱。这样"没有样本"和"低概率"能区分开，也不会出现镂空/遮挡
- `--gamma`（默认 0.6）：颜色映射提亮低概率柱，`<1` 更亮
- 灰底仅表示该格原始计数为 0；柱高仍是周期平滑后的概率

### 中心残基叠加图（20 张）

```powershell
src\process_Statistics\.venv\Scripts\python.exe src\process_Statistics\tools\plot_overlay.py `
    --probs result\static\triplet_probs.npy --out result\static\img_overlay --alpha 0.03 --floor 0.03
```

对 20 种中心残基 `m_i = X` 各出一张图，把该中心对应的 **400 个上下文**
（`m_{i-1}`, `m_{i+1}` 各 20 种）的条件分布逐层叠加：

- 每层透明度极低（`--alpha`，默认 0.03），叠加后密度自然累积出趋势
- `--floor`（默认 0.03）把低于各自峰值 3% 的格剔除，否则 400 层会在 z=0 处叠成一块不透明厚板
- `--factor 2` 时渲染降到 36×36（仅显示降采样，概率方阵不受影响）
- `--contexts N` 可只叠加 N 个随机上下文（调试用），`--only A` 只画一个中心

输出 `result\static\img_overlay\center_<X>.png`，20 张，8 进程约 41 秒。

平滑实现为可分离循环高斯核 + 2D FFT，平滑后重新归一；已验证总质量守恒（1±1e-12）
与环面平移等变性（误差 ~1e-17）。

## 实测结果（本机全量）

| 指标 | 数值 |
|---|---|
| 输入片段 | 47,794 |
| 残基 | 7,640,882 |
| 有效样本（n−2） | 7,545,294（与期望值完全一致） |
| 非空方阵 | 8000 / 8000 |
| 每三元组样本数 | 中位 714，最大 6203 |
| 单格最大概率 | 0.4545 |
| 归一化违规 | 0 |
| 统计耗时 | 8.0 s |
| 方阵文件 | 158.2 MB |
| 500 张图 | 93.6 MB，8 进程 121 s |

## 测试

```powershell
powershell -ExecutionPolicy Bypass -File src\process_Statistics\tests\run_tests.ps1
```

合成样本覆盖：分箱边界（±180° 回绕、格边界不丢点）、归一化、
稀疏矩阵保留、npy 头格式、三元组索引顺序。

## 已知限制

- 72×72 分辨率对稀有上下文偏细；本数据集 8000 个三元组全部有样本，
  但单格期望计数仍然很低，直接训练建议再做平滑或回退。
- 只读取 `segments.jsonl`，不回溯原始 PDB。
