# von_Statistics —— DPMM 源码（从 H 盘归入）

本目录是 2026-11-23 从 `H:\DATABACE\src` 归入的 **von_dpmm**（Dirichlet 过程混合 +
von Mises sine + MCMC）源码。原 `README.md` 里的路径仍是 H 盘写法，且编码有问题；
以本文件为准。

## 归入后的位置

| 原位置（H 盘） | 现位置 |
|---|---|
| `H:\DATABACE\src\{CMakeLists.txt,build.bat,config.json,src,py}` | `round\src\von_Statistics\` |
| `H:\DATABACE\init\*` | `round\datebase\date_init\`（5 个 `.bin` + `manifest.tsv` + `_done.txt`） |
| `H:\DATABACE\thr_result\von`、`von_box` | `round\result\von`、`round\result\von_box`（原本已有，哈希一致） |
| `H:\DATABACE\thr_result\staic_box`、`img`、`REPORT.md` | `round\result\staic_box`、`round\result\img`、`round\result\REPORT.md`（本次补入） |

## 配置

`config.json` 已改成相对路径（相对本目录）：

```json
"BIN_DIR" : "..\\..\\datebase\\date_init",
"MANIFEST": "..\\..\\datebase\\date_init\\manifest.tsv",
"OUT_ROOT": "..\\..\\result"
```

## 构建

`build.bat` 已从 CUDA 13.3 改为本机可用的 **CUDA 12.8 + MSVC v143**：

- `VsDevCmd.bat -vcvars_ver=14.40.33807`
- `-DCMAKE_CUDA_COMPILER=D:/Cpp/CUDA/12.8/bin/nvcc.exe`
- `-DCMAKE_CUDA_ARCHITECTURES=75`

```powershell
cd "D:\CodingPrograms\Protein Design\round\src\von_Statistics"
.\build.bat          # 增量
.\build.bat clean    # 全清重编
.\bin\von_dpmm.exe --config config.json --selftest
.\bin\von_dpmm.exe --config config.json
```

产物：`bin\von_dpmm.exe`（`bin/` 与 `build/` 已被 `.gitignore` 忽略）。

> 说明：本机 CUDA 12.8 与 VS2026 默认 v145 不兼容，所以必须用 v143；以上改动尚未实际编译验证。
> GPU 后端要求驱动 ≥580 或 CUDA 12.8 运行时；当前 573.91 驱动下程序会自动退回 CPU 后端。
