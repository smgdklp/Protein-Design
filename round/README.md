---
license: apache-2.0
pretty_name: Protein Backbone Dihedral Angles (PDB subset, phi/psi)
tags:
- protein
- pdb
- dihedral-angles
- ramachandran
- structural-biology
- protein-structure
- bioinformatics
size_categories:
- 1M<n<10M
---

# 蛋白质主链二面角数据集

## 说明

训练初期没做好程序管理,使得不能从源码编译快速重复获得数据
故把数据上传至魔塔平台

> **数据集声明**
>
> - 平台：ModelScope（公开）
> - 地址：<https://modelscope.cn/datasets/Smgdklp/protein-dihedral-dist-adjacent-aa>
> - ID：`Smgdklp/protein-dihedral-dist-adjacent-aa`
> - 中文名：**氨基酸二面角分布与相邻氨基酸关系**
> - 许可：Apache License 2.0（ModelScope 页面标注；原始 PDB 数据为 CC0）
> - 版本：v1 / 2026-11-24
>
> 本文件是代码仓库内的数据集声明；线上 README 以 ModelScope 页面为准。

从 RCSB PDB 随机抽取的 **9,349 个结构**中，提取连续、主链完整、由标准氨基酸组成的片段，
保存每个残基的 **φ / ψ 二面角 cos/sin**。面向的任务是学习

```
p(φ_i, ψ_i | m_{i-1}, m_i, m_{i+1})
```

即三肽上下文条件下的 Ramachandran 分布，可用于结构生成、序列设计打分、二面角先验等。

数据有**两套独立提取结果**，统计口径不同，**不要混用**：

| | `datebase/round_bin/`（新版） | `datebase/date_init/`（旧版 DPMM 容器） |
|---|---|---|
| 格式 | JSONL（一行一段） | 5 个 `.bin` 容器 + `manifest.tsv` |
| PDB 数 | 9,106 | 9,109 |
| 片段数 | 47,794 | 65,687 |
| 残基数 | 7,640,882 | 8,408,870 |
| 有效样本（n−2） | 7,545,294 | 8,277,496 |
| 体积 | 约 293 MB | 约 139.5 MB |

## ModelScope 仓库结构

线上仓库是“扁平”的：`round/datebase` 与 `round/result` 的内容直接放在仓库根目录，
不再保留 `datebase/`、`result/` 这两层包装。

| 线上目录 | 本地来源 | 内容 |
|---|---|---|
| `date_init/` | `round/datebase/date_init/` | 5 个 36×36 容器 `.bin` + `manifest.tsv` + `_done.txt` |
| `round_bin/` | `round/datebase/round_bin/` | `segments.jsonl`、`processed.jsonl`（9,349 个 PDB 台账） |
| `static/` | `round/result/static/` | 8000×72×72 经验直方图 `triplet_probs.npy` + 分析 txt |
| `peak_ground/` | `round/result/peak_ground/` | `ground/`（基底/残差）+ `peak/`（峰字典）+ `fin/`（拟合） |
| `final_version/` | `round/result/final_version/` | 连续核密度模型 v2：`probs.npy`、`forbidden.npy`、`peaks.jsonl` |
| `von/` | `round/result/von/` | DPMM 参数化模型（36×36） |
| `von_box/` | `round/result/von_box/` | von 参数模型的 MC 采样方阵 |
| `staic_box/` | `round/result/staic_box/` | 纯经验 36×36 方阵 |
| `README.md` | — | 线上说明（当前为 `round/result/README.md` 的内容） |
| `.gitattributes` | — | LFS 规则 |

图片类文件（`img_all/`、`imgs/`、`compare/`、`base_sheets/`、`chk_*.png`）不参与
代码仓库提交；线上仓库是否包含，以 ModelScope 文件树为准。

## 文件结构

```
datebase/
├── round_bin/
│   ├── segments.jsonl      # 293 MB，推荐入口
│   ├── processed.jsonl     # 1.8 MB，9,349 个 PDB 的提取台账
│   └── run.log
└── date_init/
    ├── 10cc8c2df62a0193.bin   # 32 MB（预分配容器）
    ├── 20ad46ce9e6d81d8.bin   # 32 MB
    ├── 47be907732671f35.bin   # 32 MB
    ├── 9e5b5a76387fc20b.bin   # 32 MB
    ├── ef3ded85438aff6a.bin   # 8.6 MB
    ├── manifest.tsv           # 2.9 MB，65,687 行片段索引
    └── _done.txt              # 9,349 个已处理 PDB ID
```

`date_init` 的 `.bin` 是**按 32 MB 预分配的容器**，一个文件里放多段；实际有效字节由
`manifest.tsv` 的 `offset` / `bytes` 指定，不能只按文件大小遍历。

## 格式 A：`round_bin/segments.jsonl`（推荐）

一行一个片段：

| 字段 | 类型 | 说明 |
|---|---|---|
| `v` | int | 格式版本，当前 1 |
| `seg` | str | 片段 id：`<PDB>_<链>_<起始残基>` |
| `pdb` | str | 来源 PDB id（大写） |
| `chain` | str | 链 id，空链为 `""` |
| `start` / `end` | int | 首末残基号 |
| `len` | int | 残基数 |
| `seq` | str | 单字母序列（长度 = `len`） |
| `cos_phi` / `sin_phi` | float[] | 与 `seq` 等长；φ 的 cos/sin |
| `cos_psi` / `sin_psi` | float[] | 与 `seq` 等长；ψ 的 cos/sin |
| `src` | str | 来源文件名 |

约定：

- φ = `C(i-1)-N-CA-C`，ψ = `N-CA-C-N(i+1)`，IUPAC 符号（α 螺旋 φ≈−60°、ψ≈−45°）；
- 端点角不存在：`cos_phi[0]=sin_phi[0]=0`、`cos_psi[n-1]=sin_psi[n-1]=0`。
  这两个 `(0,0)` 经 `atan2` 会变成假角度 0°，统计时必须跳过 `i=0` 和 `i=n-1`；
- 角度范围 `[-π, π]`，弧度。

Python 读取：

```python
import json
import numpy as np

with open("datebase/round_bin/segments.jsonl", encoding="utf-8") as f:
    for line in f:
        d = json.loads(line)
        phi = np.degrees(np.arctan2(d["sin_phi"], d["cos_phi"]))
        psi = np.degrees(np.arctan2(d["sin_psi"], d["cos_psi"]))
        # 跳过 i=0 和 i=n-1
```

`processed.jsonl` 一行一个 PDB，字段：
`pdb, status, segments, residues, dropped_segments, dropped_residues,
dist_breaks, mse, nonstandard, bytes, file, time`。
本次快照 `status=ok` 共 **9,349** 条，失败 0。

## 格式 B：`date_init/`（旧版 DPMM 容器）

`manifest.tsv`（Tab 分隔，65,687 行数据）字段：

| 列 | 说明 |
|---|---|
| `file` | 容器 `.bin` 文件名 |
| `offset` | 该片段在文件中的字节偏移 |
| `bytes` | 该片段总字节数，恒等于 `4 + 17*n` |
| `residues` | 残基数 n |
| `pdb` / `model` / `chain` | 来源信息 |

每个片段在 `offset` 处是：

```
int32 n                     残基数
n × 17 字节              每个残基一条记录
```

单条 17 字节记录（小端，无对齐）：

| 偏移 | 类型 | 含义 |
|---|---|---|
| 0 | uint8 | 氨基酸单字母 ASCII（`A`–`Y`；非标准为 `X`） |
| 1 | float32 | cos φ |
| 5 | float32 | sin φ |
| 9 | float32 | cos ψ |
| 13 | float32 | sin ψ |

缺失角用 **NaN** 哨兵：段首 φ 为 `(NaN, NaN)`，段尾 ψ 为 `(NaN, NaN)`。

Python 读取（按 manifest 逐段）：

```python
import struct
import numpy as np

with open("datebase/date_init/manifest.tsv", encoding="utf-8") as f:
    header = f.readline()
    for line in f:
        file, offset, nbytes, n, pdb, model, chain = line.rstrip("\n").split("\t")
        with open(f"datebase/date_init/{file}", "rb") as bf:
            bf.seek(int(offset))
            n = struct.unpack("<i", bf.read(4))[0]
            raw = bf.read(17 * n)
        vals = np.array(
            [struct.unpack_from("<4f", raw, 17 * i + 1) for i in range(n)],
            dtype=np.float32,
        )
        aa = np.frombuffer(raw, dtype=np.uint8).reshape(n, 17)[:, 0]
        phi = np.degrees(np.arctan2(vals[:, 1], vals[:, 0]))
        psi = np.degrees(np.arctan2(vals[:, 3], vals[:, 2]))
```

## 提取规则

两套数据都来自同一批 9,349 个 PDB，但由不同版本的提取器生成，因此片段数/残基数不同：

- 只读第一个 `MODEL`（NMR 多模型不重复计数）；
- 20 种标准氨基酸；`MSE` 视作 `M`；水、配体、核酸忽略；
- 在 `TER`、换链、缺 N/CA/C 任一主链原子、相邻残基 `C–N` 距离 > 2.5 Å 处切段；
- 段长 < 3 丢弃；
- 二面角符号与 Biopython / PyMOL 一致，不做翻转。

`date_init` 是较早的容器格式，保留了更多片段（65,687 段）；`round_bin` 是较严格的
标准氨基酸口径（47,794 段）。需要单一口径时选 `round_bin`，不要把两者拼接。

## 上下文索引（下游常用）

- 字母表：`ACDEFGHIKLMNPQRSTVWY`
- 三元组：`m_{i-1}, m_i, m_{i+1}`，共 `20³ = 8000` 个上下文
- `idx = (a*20 + b)*20 + c`，`a/b/c` 为 0..19
- 常用分箱：72×72（5°/格）或 36×36（10°/格）

## 已知偏置与限制

- **有序化偏置**：只保留主链完整、无断链的残基，无序区/柔性区系统性缺席；
- **同源冗余**：7 万多段来自 9 千个 PDB，同源结构/同一蛋白多段都在内。做显著性检验或
  训练/测试切分时必须按 **PDB id 或链/家族分块**，不能按残基随机切；
- 稀疏上下文：C / W 为中心的某些三元组样本很少，直接拟合会过拟合；
- 两套数据口径不同（见上），不要混用；
- 本数据集不提供官方 train/test split。

## 许可与来源

- 原始结构来自 [RCSB PDB](https://www.rcsb.org/)，数据许可 **CC0**；
- 本数据集在 ModelScope 页面上标注为 **Apache License 2.0**（派生提取部分按此发布）；
- 不提供任何正确性担保，使用前请自行核对几何与符号约定。

## 版本

- `v1` / 2026-11-24：从 9,349 个 PDB 的本地快照导出；两套格式同时提供。
- 原始 PDB 清单见 `datebase/round_bin/processed.jsonl` 与 `datebase/date_init/_done.txt`。
