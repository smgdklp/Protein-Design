#pragma once
// .pdb / .pdb.gz -> 定长二面角二进制
//
// 落盘结构 struct List(小端, 紧凑无填充):
//   [int32 _len]                    残基数 N
//   ( [char _kind][float32 ca][float32 sa][float32 cb][float32 sb] ) x N   // 每单元 17 字节
//   文件总长 = 4 + 17 * N
//
// ca/sa = cos φ / sin φ, cb/sb = cos ψ / sin ψ;
// 不存在或算不出的角(序列两端、链断裂、骨架缺失)四个分量一律写 0.0f。
#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

namespace pbd {

// _kind: 1..20 = 标准氨基酸(三字缩写字母序), 21 = 非标准 / 缺 N-CA-C 骨架
inline constexpr char kKindUnknown = 21;

inline constexpr std::size_t kUnitBytes   = 17; // char + 4 x float32
inline constexpr std::size_t kHeaderBytes = 4;  // int32 _len

struct Residue {
    char kind = kKindUnknown;
    unsigned char mask = 0;  // bit0 N, bit1 CA, bit2 C
    unsigned char brk = 0;   // 1 = 与前一残基不连续(换链 / TER / 编号跳号)
    char chain = ' ';
    std::int32_t seq = 0;
    char icode = ' ';
    float n[3] = {0.0f, 0.0f, 0.0f};
    float ca[3] = {0.0f, 0.0f, 0.0f};
    float c[3] = {0.0f, 0.0f, 0.0f};
};

// 三字缩写 -> 1..20; 认不出(非标准/修饰残基/其它)一律 21
char KindOf(std::string_view three) noexcept;

// 后缀判断(ASCII 大小写不敏感); suffix 必须是全小写
bool EndsWithNoCase(std::string_view s, std::string_view lower_suffix) noexcept;

// abc.pdb / abc.pdb.gz -> abc.bin(只换后缀, 文件名其余部分保留)
std::filesystem::path BinNameFor(const std::filesystem::path& source);

// ---------------------------------------------------------------- 内存池: PDB 侧
// 只负责把文本读进来并折成残基数组; 缓冲区跨文件复用, 不反复分配。
class pbd_pool {
public:
    bool Load(const std::filesystem::path& file, std::string& err);
    bool Parse(std::string& err);

    const std::vector<Residue>& residues() const noexcept { return res_; }
    std::size_t text_bytes() const noexcept { return text_.size(); }
    std::size_t capacity_bytes() const noexcept;

    void Clear();

private:
    std::vector<char> text_;   // 文件正文 / 解压后正文
    std::vector<Residue> res_; // 残基序列
};

// ---------------------------------------------------------------- 内存池: 二进制侧
// 只攒字节, 不碰文件系统细节之外的逻辑。
class bin_pool {
public:
    void Reset();
    void Append(char kind, float cos_phi, float sin_phi, float cos_psi, float sin_psi);
    bool Save(const std::filesystem::path& out, std::string& err) const;

    std::size_t residues() const noexcept { return residues_; }
    std::size_t payload_bytes() const noexcept { return bytes_.size() - kHeaderBytes; }
    std::size_t capacity_bytes() const noexcept { return bytes_.capacity(); }

private:
    void AppendF32(float v);

    std::vector<unsigned char> bytes_; // 前 kHeaderBytes 字节留给 _len
    std::size_t residues_ = 0;
};

// ---------------------------------------------------------------- 单文件处理单元
// 一个 maker 一个线程, 内部两个池都是私有的, 全程无锁。
class SimgleAnle_maker {
public:
    void Init(std::string bin_path); // 只接 BIN_PATH

    // 只接单个 .pdb / .pdb.gz 路径
    bool Work(const std::filesystem::path& file,
              std::uint32_t& residue_count,
              std::uint64_t& out_bytes,
              std::string& err);

    std::size_t capacity_bytes() const noexcept;

private:
    std::filesystem::path bin_dir_;
    pbd_pool pbd_;
    bin_pool bin_;
};

} // namespace pbd
