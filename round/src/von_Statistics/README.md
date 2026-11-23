# von_dpmm — 二面角条件分布 P(φ,ψ | m₋₁,m₀,m₊₁)

从 `H:\DATABACE\init` 的 `cpp_round` 二进制里,对 20³ = 8000 个三元组上下文各拟合一个
**Dirichlet 过程混合(von Mises sine 成分)**,按 Lennox(2009)+`高斯.md` 的先验/迭代方案跑 MCMC,
把后验预测密度落在 36×36 网格上,再压缩回参数;最后从参数化模型做离散采样得到方阵。

## 构建

CUDA 13.3 没有注册到 VS 的 BuildCustomizations,所以**不能用 "Visual Studio" 生成器**;
`build.bat` 用 NMake Makefiles(CMake 直接调 nvcc),主机编译器仍是 MSVC 14.51。

```powershell
cd H:\DATABACE\src
.\build.bat          # 增量
.\build.bat clean    # 全清重编
```

产物:`bin\von_dpmm.exe`。

## 运行

```powershell
cd H:\DATABACE\src
.\bin\von_dpmm.exe --config config.json          # 全流程
.\bin\von_dpmm.exe --config config.json --selftest   # 数学自检
.\bin\von_dpmm.exe --only-hist                   # 只做经验直方图
.\bin\von_dpmm.exe --no-sample                   # 只到参数
```

## 后端

程序先探测 GPU;探测失败自动退回 **CPU 后端**(与 kernel 同算法、同随机数流,多线程跑上下文)。

- 本机实测:`Quadro RTX 3000 / driver 573.91` 只支持到 **CUDA 12.8**,而唯一装着的工具链是
  **CUDA 13.3** → `cudaGetDevice` 报错,GPU 后端不可用。修法二选一:
  1. 升级 NVIDIA 驱动到 ≥ 580(M 系列之后),CUDA 13.3 即可用;
  2. 装 CUDA 12.8 toolkit,改 `CMakeLists.txt`/`build.bat` 里的 nvcc 路径,并把
     `CMAKE_CUDA_ARCHITECTURES` 保持 75。
  两者都要下载,按本机规矩必须请示。
- CPU 后端跑全量(8000 上下文、Tier A 2 链 6000 迭代)约 80–100 分钟;GPU 后端设计上是分钟级。

## 配置项(config.json)

| 键 | 含义 |
|---|---|
| `BIN_DIR` | `cpp_round` 输出的 `*.bin` 目录 |
| `OUT_ROOT` | 四个产物目录的父目录 |
| `SEED` | 全部随机流的种子(改了就换一套 MCMC 样本) |
| `TIER_A_N` / `TIER_B_N` | 分级阈值(默认 2000 / 200) |
| `A_ITERS/A_BURN/A_THIN/A_CHAINS` | Tier A 预算(默认 6000/1000/10/2) |
| `B_...` / `C_...` | Tier B / C 预算 |
| `TAU0` | DP 质量参数(默认 1.0) |
| `K10/K20/LAM0` | H1 = vM sine 先验的精度(默认 1/π²) |
| `WISHART_NU0/WISHART_B0` | H2 = Wishart 先验(默认 2 / 0.25) |
| `K_AUX` | Algorithm 8 每轮辅助簇数(默认 3) |
| `K_MAX_COMP` | 压缩后最大成分数(默认 8) |
| `M_SAMPLES` | von_box 每上下文采样数(默认 200000) |
| `MAX_CTX` | 只跑前 N 个上下文(调试) |

## 结构

```
src/main.cpp    编排、分级统计、一致性核对、--selftest
src/model.h     常量/索引/二进制格式/接口
src/config.cpp  容错 JSON
src/io.cpp      bin 扫描、CSR、Bessel 表、二进制写出
src/cpu.cpp     CPU 后端(参考实现):MCMC、式(7) 网格、压缩、采样
src/gpu.cu      CUDA 核:直方图、MCMC、压缩、采样
src/py/plot.py  画图(随机 5 个上下文 × 三组,环面插值)
```

## 算法要点(与论文的对应)

- 成分密度(论文式 1):`p(φ,ψ|μ,ν,Ω) = C·exp{κ₁cos(φ−μ) + κ₂cos(ψ−ν) + λ sin(φ−μ)sin(ψ−ν)}`,
  `Ω = [[κ₁,−λ],[−λ,κ₂]]`;`λ² < κ₁κ₂` 等价于 Ω 正定,即单峰。归一化常数 C 用 Bessel 级数,
  GPU/CPU 都走"log(e^{-κ}I_m(κ)) 查表 + logsumexp"。
- H1 = vM sine(μ₀=ν₀=0, κ₁₀=κ₂₀=1/π², λ₀=0);H2 = Wishart₂(ν₀=2, B₀=0.25·I);τ₀=1。
- 分配步:Neal(2000)Algorithm 8 的**块更新变体**(一轮 sweep 内辅助簇固定,点并行独立分配),
  等价于对 DP 的截断近似;空簇每轮回收。
- 均值步:`μ|ν`、`ν|μ` 的全条件**恰好是单变量 von Mises**,所以用精确两步 Gibbs
  (论文用二维拒绝采样,这里换成等价但更适合并行的形式)。
- 精度步:`(κ₁,κ₂,λ)` 用随机游走 MH,先验取 Wishart₂(论文的 Wishart 直采是伪共轭近似,
  这里用 MH 保证采的是正确后验)。
- 密度估计(论文式 7):`p̂ = (1/B)Σ_b p(·|μ^b,ν^b,Ω^b)`,在 36×36 格心求值 × 格面积,再归一。
- 压缩:用 EM 把 `p̂` 拟合成 K'≤8 个 vM sine 成分,K' 由 BIC 选;初值取 MCMC 末态簇。

## 验收

`--selftest` 会做三件事,任一不过就不该信后面的数字:

1. Bessel 表插值 vs 直算级数(log 尺度);
2. `λ=0` 时 `log C` 必须等于 `−log(4π²I₀(κ₁)I₀(κ₂))`;
3. 对 κ₁=20,κ₂=30,λ=8 在 720×720 上数值积分,密度积分必须 = 1。

实测:插值最大偏差 `5.9e-06`,`log C` 偏差 `0`,`数值积分 1.000000`。
