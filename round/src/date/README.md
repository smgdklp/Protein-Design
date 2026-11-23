# pdbextract —— PDB 标准氨基酸片段提取器

从 PDB 文件中提取**连续、主链完整、且全部由 20 种标准氨基酸组成**的片段，
按残基输出序列、单字母标签以及 φ/ψ 二面角的 cos/sin。

## 构建

```powershell
cmake -S "D:\CodingPrograms\Protein Design\round" -B "D:\CodingPrograms\Protein Design\round\build" -G "Visual Studio 18 2026" -A x64
cmake --build "D:\CodingPrograms\Protein Design\round\build" --config Release --parallel
```

产物：`round\bin\Release\pdbextract.exe`，`pbd_config.json` 会在构建后自动复制到 exe 同目录。

## 配置（`pbd_config.json`，与 exe 同目录）

```json
{
  "PBD_PATH": "D:\\Huggingface\\databace\\random_pbd",
  "BIN_PATH": "D:\\CodingPrograms\\Protein Design\\round\\datebase\\round_bin",
  "NUM": 5
}
```

| 键 | 含义 | 缺省 |
|---|---|---|
| `PBD_PATH`（也接受 `PDB_PATH`） | 输入 PDB 文件夹 | 必填 |
| `BIN_PATH` | 输出文件夹，不存在会自动创建 | 必填 |
| `NUM` | 线程数（数字或数字字符串均可，限制 1–64） | 5 |
| `BREAK_DIST` | 相邻残基 C(i-1)–N(i) 距离阈值（Å），超过即切分；`0` 关闭 | 2.5 |
| `MIN_LEN` | 片段最小残基数，不足整段丢弃 | 3 |
| `RECURSIVE` | 是否递归扫描子目录 | false |

解析器是**宽松模式**：允许尾逗号、缺逗号、`//` 与 `/* */` 注释；
字符串里只有 `\"`、`\\`、`\/`、`\uXXXX` 会被转义，其它反斜杠原样保留，
所以 `"D:\round"` 这种没转义的 Windows 路径也能直接用。

## 提取规则

1. 只读第一个 `MODEL`（无 MODEL 记录则读全文）。
2. 残基白名单 = 20 种标准氨基酸；`MSE`（硒代甲硫氨酸，通常记为 HETATM）**保留**并按 `MET` 处理，
   序列字符为 `M`。其它 HETATM（水、配体）忽略，DAT/RNA（ATOM 记录）视为非标准。
3. 遇到下列情况切分片段：非标准残基、`TER`、链切换、缺 N/CA/C 任一主链原子、
   相邻残基 `C–N` 距离 > `BREAK_DIST`。
4. 片段长度 < `MIN_LEN` 时整段丢弃（台账中记入 `dropped_segments`）。
5. 二面角只算 **φ（C(i-1)-N-CA-C）** 与 **ψ（N-CA-C-N(i+1)）**，存 cos/sin。
   符号约定与 Biopython / PyMOL 一致（α 螺旋 φ≈-60°、ψ≈-45°）。
   约定 `cos_phi[0]=sin_phi[0]=0`、`cos_psi[n-1]=sin_psi[n-1]=0`（端点角不存在），
   其余位置保证有效；`MIN_LEN=3` 恰好保证每段至少一个残基 φ、ψ 都有定义。

altLoc 处理：同名原子按优先级 `空 > A > 其它` 取，优先级相同取 occupancy 大者。

## 输出

### `segments.jsonl`（一行 = 一个片段，同时充当交叉索引）

```json
{"v":1,"seg":"1ABC_A_12","pdb":"1ABC","chain":"A","start":12,"end":290,"len":279,
 "seq":"MSTN...","cos_phi":[...],"sin_phi":[...],"cos_psi":[...],"sin_psi":[...],
 "src":"1ABC.pdb"}
```

| 字段 | 说明 |
|---|---|
| `seg` | 片段唯一 id：`<PDB>_<链>_<起始残基号>[插入码]`，空链用 `_` |
| `pdb` | 来源 PDB id（文件名去扩展名、大写） |
| `chain` | 链 ID，空链为 `""` |
| `start` / `end` | 片段首末残基号 |
| `seq` | 单字母序列（即逐残基标签，char） |
| `cos_phi` … `sin_psi` | 与 `seq` 等长的稠密数组，6 位小数 |
| `src` | 来源文件名 |

### `processed.jsonl`（台账，断点续跑依据）

```json
{"pdb":"1ABC","status":"ok","segments":7,"residues":1832,"dropped_segments":1,
 "dropped_residues":2,"dist_breaks":0,"mse":1,"nonstandard":4,"bytes":1048576,
 "file":"1ABC.pdb","time":"2026-11-08T13:45:02"}
```

失败文件记 `status:"failed"` 并带 `error`，下次运行会重试。

### `run.log`

运行日志：配置、文件数、进度、汇总。

## 用法

```
pdbextract [--config <pbd_config.json>] [--limit N] [--threads N] [--rebuild] [--stats]
```

- `--config` 缺省依次找 exe 同目录、当前目录的 `pbd_config.json`
- `--limit N` 只处理前 N 个待办文件（冒烟用）
- `--rebuild` 删除已有 `segments.jsonl` / `processed.jsonl` 后全量重跑
- `--stats` 统计已有产出的片段数、残基数、角度自洽性、φ/ψ 直方图

**幂等**：`processed.jsonl` 中 `status=ok` 的 id 会被跳过；重复运行不会重复产出。

## 测试

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\src\date\tests\run_tests.ps1
```

脚本会在 `build\test_runs\` 下生成合成 PDB（标准链、短片段、非标准残基、MSE、
几何断链、缺主链原子、多模型），并断言切分结果、序列、长度、角度自洽与宽松 JSON 解析。

## 已知限制

- 只认 PDB 文本格式的 `ATOM`/`HETATM`，不支持 mmCIF、gzip 压缩包。
- 坐标按 PDB 固定列宽解析；超过 9999 残基的 hybrid-36 编号只会影响 `start/end` 的显示，
  不影响切分与序列。
- 未做几何质量过滤（键长/键角异常），只做距离断链。
