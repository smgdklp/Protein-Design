# cpp_pbddatemaker

把 `.pdb` / `.pdb.gz` 里的氨基酸序列与主链二面角，压成定长二进制 `.bin`。
**本文件是内存结构、`_kind` 编码、哈希契约的唯一权威声明**，改这三样必须同步改这里。

## 一、构建与运行

```powershell
# 配置(生成器固定 Visual Studio 18 2026 / x64, 不用 Ninja)
cmake -S "D:\CodingPrograms\Protein Design\tool\src\cpp_pbddatemaker" `
      -B "D:\CodingPrograms\Protein Design\tool\src\cpp_pbddatemaker\build" `
      -G "Visual Studio 18 2026" -A x64

# 编译
cmake --build "D:\CodingPrograms\Protein Design\tool\src\cpp_pbddatemaker\build" --config Release

# 运行(配置从 exe 同目录读, 不看当前工作目录)
"D:\CodingPrograms\Protein Design\tool\result\pbd_datamaker\pbddatemaker.exe"
```

构建时 `config.json` 被复制到 `tool/result/pbd_datamaker/`（与 exe 同目录）。

## 二、config.json（全部参数都是字符串）

```json
{
    "PBD_PATH" : "D:\\Huggingface\\databace\\random_pbd",
    "BIN_PATH" : "D:\\Huggingface\\databace\\random_pbd\\bin",
    "NUM" : "5",
    "OUTPUT" : "SimgleAnle"
}
```

| 键 | 含义 |
|---|---|
| `PBD_PATH` | 输入根目录，递归遍历，只收 `.pdb` / `.pdb.gz`；`BIN_PATH` 子树会被剪掉 |
| `BIN_PATH` | 输出目录，自动创建；文件名保留主名，只换后缀 `10GN.pdb` → `10GN.bin` |
| `NUM` | worker 数；每个 worker 一个 `SimgleAnle_maker` + 一个独占线程 |
| `OUTPUT` | 分流依据，目前只认 `"SimgleAnle"`；不认识直接报错退出 |

### 哈希契约（跨语言必须一致）

所有键、所有枚举取值对外都是字符串，进入程序后立即折成 `uint64`，之后内部只比较哈希码：

```
算法  FNV-1a 64
输入  UTF-8 字节, 区分大小写, 不做归一化/去空白
offset 0xCBF29CE484222325
prime  0x100000001B3
```

与 `tool/src/cpp_pdbtrans` 使用同一算法，同一字符串在两边得到同一哈希码。
`NUM` 除了哈希身份之外，还会用十进制解析成真实数字（哈希不参与算术）。

## 三、磁盘结构 `struct List`

小端，紧凑无填充，一个输入文件对应一个 `.bin`：

```
偏移 0        int32   _len                 残基数 N
偏移 4+17k    char    _kind                第 k 个残基 (k = 0..N-1)
偏移 5+17k    float32 ca = cos φ
偏移 9+17k    float32 sa = sin φ
偏移 13+17k   float32 cb = cos ψ
偏移 17+17k   float32 sb = sin ψ
文件总长 = 4 + 17 * N
```

即 `char` 与 4 个 `float32` 交替，每单元 **17 字节**。17 不是 4 的倍数，所以单元内的 float
**不对齐**：读写一律按字节 `memcpy`，不要取地址当 `float*` 用。

### `_kind` 编码（国际标准顺序：三字缩写字母序）

| `_kind` | 三字 | 单字 | `_kind` | 三字 | 单字 |
|---|---|---|---|---|---|
| 1 | ALA | A | 11 | LEU | L |
| 2 | ARG | R | 12 | LYS | K |
| 3 | ASN | N | 13 | MET | M |
| 4 | ASP | D | 14 | PHE | F |
| 5 | CYS | C | 15 | PRO | P |
| 6 | GLN | Q | 16 | SER | S |
| 7 | GLU | E | 17 | THR | T |
| 8 | GLY | G | 18 | TRP | W |
| 9 | HIS | H | 19 | TYR | Y |
| 10 | ILE | I | 20 | VAL | V |
| — | — | — | **21** | **非标准/统一未知** | **X** |

`21` 同时覆盖两种情况：非 20 种标准残基（MSE、SEC、UNK、各种修饰残基），以及骨架原子缺失
（有 Cα 但缺 N 或 C）的残基。**这类残基不会被丢弃，仍占一个单元、计入 `_len`。**
`_kind` 只表示分类，不是几何开关：骨架完整的非标准残基照样算出真实 φ/ψ 写入。

## 四、残基筛选与角语义

- 只吃第一个 `MODEL`；`ENDMDL` 即停。
- 一个残基进入序列的条件：**有 Cα，且至少有 N 或 C**。
  这条规则天然排除水、离子、单原子配体，以及 DNA/RNA（它们没有 Cα，不是氨基酸残基）。
- altloc 优先级：空白 > `A` > `B` > `C` > 其它（不按占有率挑；同一残基的 N/CA/C 可能来自同一 altloc）。
- 链断裂判定：换链、遇到 `TER`、残基编号跳号（同号但插入码更大视为连续）→ 标记 `brk`。
- `φ = dihedral(C[i-1], N[i], CA[i], C[i])`，`ψ = dihedral(N[i], CA[i], C[i], N[i+1])`，
  符号为 IUPAC 标准（已用 PDB 自带 `REMARK 500` 官方 φ/ψ 逐个对过，见第七节）。
- **不存在的角一律写字面 `0.0f`**（不是 `(cos,sin)=(1,0)`，也不是 NaN），两个角各自独立判定：
  - φ 不存在 → `ca = sa = 0`：序列首端、`brk` 之前、本残基缺 N/CA/C、前一残基缺 C。
  - ψ 不存在 → `cb = sb = 0`：序列末端、`brk` 之后、本残基缺 N/CA/C、后一残基缺 N。
  - 读端判据就是「该角的 cos 与 sin 同时为 0」；两者不会同时为 0 于合法角（cos²+sin²=1）。

## 五、读取示例

```cpp
#include <cstdint>
#include <cstring>
#include <fstream>
#include <vector>

struct Unit { char kind; float ca, sa, cb, sb; }; // 17 字节, 但别直接 cast 内存

std::vector<Unit> ReadBin(const char* path) {
    std::ifstream f(path, std::ios::binary);
    std::int32_t n = 0;
    f.read(reinterpret_cast<char*>(&n), 4);

    std::vector<Unit> out(static_cast<std::size_t>(n));
    unsigned char raw[17];
    for (std::int32_t i = 0; i < n; ++i) {
        f.read(reinterpret_cast<char*>(raw), 17);
        Unit u{};
        u.kind = static_cast<char>(raw[0]);
        std::memcpy(&u.ca, raw + 1, 4);
        std::memcpy(&u.sa, raw + 5, 4);
        std::memcpy(&u.cb, raw + 9, 4);
        std::memcpy(&u.sb, raw + 13, 4);
        out[static_cast<std::size_t>(i)] = u;
    }
    return out;
}
```

## 六、目录

```
CMakeLists.txt          构建定义(产物落在 tool/result/pbd_datamaker)
config.json             唯一配置源
src/main.cpp            Mian: 读配置 / OUTPUT 哈希分流 / NUM 线程 / UTF-8 进度
src/config.hpp|cpp      配置读取(扁平 JSON -> 哈希表)
src/hash.hpp            FNV-1a 64 统一哈希
src/json_lite.hpp       极简 JSON(扁平对象)
src/utf8path.hpp        UTF-8 <-> filesystem::path
src/reader.hpp|cpp      pbd_pool / bin_pool / SimgleAnle_maker
```

## 七、实测与校验（2026 全量跑通）

按仓库内 `config.json` 对 `D:\Huggingface\databace\random_pbd` 全量运行：

| 项 | 结果 |
|---|---|
| 输入 | 4679 个 `.pdb`（无 `.gz`），5 线程，用时 ≈ 4.2 s |
| 输出 | 4679 个 `.bin`，共 64,276,710 字节 |
| 残基总数 | 3,779,882（`kind = 1..20` 3,770,271 + `kind = 21` 9,611） |
| 空结构 | 106 个纯 DNA/RNA 条目 → 0 残基，落 4 字节文件（`_len = 0`） |
| 长度自洽 | 4679/4679 全部满足 `文件大小 == 4 + 17 × _len`，失败 0 |
| 失败 | 0 |

**真值校验**：用 PDB 自带的 `REMARK 500`（官方 Kleywegt 校验算出的 φ/ψ 异常表）反查
前 40 个文件里 10 个含该段的条目，共比 290 个角：

- 288 个吻合到 0.05° 量级（PDB 表格只有两位小数）。
- 2 个偏差来自 altloc 选择：`10PW` 的 `E174`（PDB 用的是 altloc **B**，本工程按文档规则取 **A**）
  与 `F84`（链 `F` 的 `LYS83` 同为 altloc 分歧，差 0.9°）。属策略差异，不是算错。

> 注意：兄弟工程 `tool/src/cpp_pdbtrans` 旧产物里的 `sa`/`sb` 是反号的（它自己把 atan2 的
> 分子取反了），cos 分量与本工程一致。两者**不能**直接混用，也不要拿它当 φ/ψ 真值。
