#include "reader.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstring>
#include <fstream>

#include "utf8path.hpp"
#include "zlib.h"

namespace pbd {

namespace {

// PDB 定长列之外的部分一律按空白裁掉
std::string_view Field(const char* line, std::size_t len, std::size_t off, std::size_t width) noexcept {
    if (off >= len) return {};
    std::size_t b = off;
    std::size_t e = (off + width < len) ? off + width : len;
    while (b < e && (line[b] == ' ' || line[b] == '\t')) ++b;
    while (e > b && (line[e - 1] == ' ' || line[e - 1] == '\t' || line[e - 1] == '\r')) --e;
    return std::string_view(line + b, e - b);
}

bool ParseF(std::string_view s, float& out) noexcept {
    if (s.empty()) return false;
    const auto r = std::from_chars(s.data(), s.data() + s.size(), out);
    return r.ec == std::errc{} && r.ptr == s.data() + s.size();
}

std::int32_t ParseI(std::string_view s) noexcept {
    if (s.empty()) return 0;
    std::int32_t v = 0;
    const auto r = std::from_chars(s.data(), s.data() + s.size(), v);
    return r.ec == std::errc{} ? v : 0;
}

// altLoc 优先级: 空白 > A > B > C > 其它
unsigned char AltRank(char alt) noexcept {
    switch (alt) {
    case ' ': return 4;
    case 'A': return 3;
    case 'B': return 2;
    case 'C': return 1;
    default:  return 0;
    }
}

// 二面角(弧度), 标准 IUPAC 符号:p0-p1-p2-p3 绕 p1->p2 轴
float Dihedral(const float p0[3], const float p1[3], const float p2[3], const float p3[3]) noexcept {
    const double b0[3] = {static_cast<double>(p0[0]) - p1[0],
                          static_cast<double>(p0[1]) - p1[1],
                          static_cast<double>(p0[2]) - p1[2]};
    double b1[3] = {static_cast<double>(p2[0]) - p1[0],
                    static_cast<double>(p2[1]) - p1[1],
                    static_cast<double>(p2[2]) - p1[2]};
    const double b2[3] = {static_cast<double>(p3[0]) - p2[0],
                          static_cast<double>(p3[1]) - p2[1],
                          static_cast<double>(p3[2]) - p2[2]};

    const double bl = std::sqrt(b1[0] * b1[0] + b1[1] * b1[1] + b1[2] * b1[2]);
    if (!(bl > 1e-9)) return 0.0f;
    for (double& v : b1) v /= bl;

    const double d0 = b0[0] * b1[0] + b0[1] * b1[1] + b0[2] * b1[2];
    const double d2 = b2[0] * b1[0] + b2[1] * b1[1] + b2[2] * b1[2];

    const double v[3] = {b0[0] - d0 * b1[0], b0[1] - d0 * b1[1], b0[2] - d0 * b1[2]};
    const double w[3] = {b2[0] - d2 * b1[0], b2[1] - d2 * b1[1], b2[2] - d2 * b1[2]};

    const double x = v[0] * w[0] + v[1] * w[1] + v[2] * w[2];
    const double cx = b1[1] * v[2] - b1[2] * v[1];
    const double cy = b1[2] * v[0] - b1[0] * v[2];
    const double cz = b1[0] * v[1] - b1[1] * v[0];
    const double y = cx * w[0] + cy * w[1] + cz * w[2];

    if (x == 0.0 && y == 0.0) return 0.0f;
    return static_cast<float>(std::atan2(y, x));
}

} // namespace

bool EndsWithNoCase(std::string_view s, std::string_view lower_suffix) noexcept {
    if (s.size() < lower_suffix.size()) return false;
    s.remove_prefix(s.size() - lower_suffix.size());
    for (std::size_t i = 0; i < lower_suffix.size(); ++i) {
        char a = s[i];
        if (a >= 'A' && a <= 'Z') a = static_cast<char>(a - 'A' + 'a');
        if (a != lower_suffix[i]) return false;
    }
    return true;
}

char KindOf(std::string_view t) noexcept {
    if (t.size() != 3) return kKindUnknown;
    // 国际标准顺序: 三字缩写字母序 ALA..VAL -> 1..20
    static constexpr char kTable[20][4] = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"};
    for (int i = 0; i < 20; ++i) {
        if (t[0] == kTable[i][0] && t[1] == kTable[i][1] && t[2] == kTable[i][2]) {
            return static_cast<char>(i + 1);
        }
    }
    return kKindUnknown;
}

std::filesystem::path BinNameFor(const std::filesystem::path& source) {
    std::string name = Utf8FromPath(source.filename());
    if (EndsWithNoCase(name, ".gz")) name.resize(name.size() - 3u);
    if (EndsWithNoCase(name, ".pdb")) name.resize(name.size() - 4u);
    name += ".bin";
    return PathFromUtf8(name);
}

// ---------------------------------------------------------------- pbd_pool

bool pbd_pool::Load(const std::filesystem::path& file, std::string& err) {
    text_.clear();
    err.clear();

    const std::string name = Utf8FromPath(file.filename());
    if (EndsWithNoCase(name, ".gz")) {
        gzFile f = gzopen_w(file.wstring().c_str(), "rb");
        if (f == nullptr) {
            err = "gzopen 失败: " + Utf8FromPath(file);
            return false;
        }
        gzbuffer(f, 1u << 18);
        char buf[1 << 16];
        for (;;) {
            const int n = gzread(f, buf, static_cast<unsigned>(sizeof buf));
            if (n < 0) {
                int zerr = 0;
                const char* msg = gzerror(f, &zerr);
                err = std::string("gzread 失败: ") + (msg != nullptr ? msg : "") +
                      " (zlib err=" + std::to_string(zerr) + ")";
                gzclose(f);
                return false;
            }
            if (n == 0) break;
            text_.insert(text_.end(), buf, buf + n);
        }
        int zerr = 0;
        gzerror(f, &zerr);
        const bool bad = (zerr != Z_OK && zerr != Z_STREAM_END);
        gzclose(f);
        if (bad) {
            err = "gz 解压不完整";
            return false;
        }
        return true;
    }

    std::ifstream in(file, std::ios::binary);
    if (!in) {
        err = "打不开输入文件: " + Utf8FromPath(file);
        return false;
    }
    in.seekg(0, std::ios::end);
    const std::streampos pos = in.tellg();
    if (pos > 0) {
        text_.resize(static_cast<std::size_t>(pos));
        in.seekg(0, std::ios::beg);
        in.read(text_.data(), static_cast<std::streamsize>(pos));
        if (!in) {
            err = "读输入文件失败: " + Utf8FromPath(file);
            return false;
        }
    }
    return true;
}

bool pbd_pool::Parse(std::string& err) {
    res_.clear();
    err.clear();

    struct Pending {
        bool active = false;
        bool ter = false;
        char name[3] = {' ', ' ', ' '};
        std::int32_t seq = 0;
        char chain = ' ';
        char icode = ' ';
        unsigned char mask = 0;      // bit0 N, bit1 CA, bit2 C
        unsigned char rank[3] = {0, 0, 0};
        float xyz[3][3] = {{0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}};
    };

    // 只有"带 Cα 且至少还有一个骨架原子"的残基才算氨基酸残基。
    // 这条规则天然排除水、离子、配体、以及 DNA/RNA(P/O5'/C4' 没有 Cα)。
    const auto commit = [this](Pending& cur) {
        if (!cur.active) return;
        cur.active = false;
        const bool has_ca = (cur.mask & 0x02u) != 0;
        const bool has_anchor = (cur.mask & 0x05u) != 0; // N 或 C
        if (!has_ca || !has_anchor) return;

        Residue r;
        r.kind = KindOf(std::string_view(cur.name, 3));
        r.mask = cur.mask;
        r.chain = cur.chain;
        r.seq = cur.seq;
        r.icode = cur.icode;
        for (int i = 0; i < 3; ++i) {
            r.n[i] = cur.xyz[0][i];
            r.ca[i] = cur.xyz[1][i];
            r.c[i] = cur.xyz[2][i];
        }
        r.brk = cur.ter ? 1u : 0u;
        if (!res_.empty()) {
            const Residue& prev = res_.back();
            if (prev.chain != r.chain) {
                r.brk = 1;
            } else if (r.seq != prev.seq + 1 &&
                       !(r.seq == prev.seq && r.icode > prev.icode)) {
                r.brk = 1;
            }
        }
        res_.push_back(r);
    };

    Pending cur;
    bool pending_ter = false;
    int model = 0;

    const char* p = text_.data();
    const char* const end = p + text_.size();
    while (p < end) {
        const char* nl = static_cast<const char*>(std::memchr(p, '\n', static_cast<std::size_t>(end - p)));
        const char* const eol = (nl != nullptr) ? nl : end;
        const std::size_t len = static_cast<std::size_t>(eol - p);

        const std::string_view rec = Field(p, len, 0, 6);
        if (rec == "MODEL") {
            ++model;
            if (model > 1) break; // 只吃第一个 model
        } else if (rec == "ENDMDL") {
            break;
        } else if (rec == "TER") {
            commit(cur);
            pending_ter = true;
        } else if (rec == "ATOM" || rec == "HETATM") {
            const std::string_view res_name = Field(p, len, 17, 3);
            char rn[3] = {' ', ' ', ' '};
            const std::size_t rn_len = (res_name.size() < 3u) ? res_name.size() : 3u;
            if (rn_len > 0) std::memcpy(rn, res_name.data(), rn_len);

            const std::int32_t seq = ParseI(Field(p, len, 22, 4));
            const char chain = (len > 21u) ? p[21] : ' ';
            const char icode = (len > 26u) ? p[26] : ' ';
            const char alt = (len > 16u) ? p[16] : ' ';

            const bool same = cur.active && cur.seq == seq && cur.chain == chain &&
                              cur.icode == icode && cur.name[0] == rn[0] &&
                              cur.name[1] == rn[1] && cur.name[2] == rn[2];
            if (!same) {
                commit(cur);
                cur = Pending{};
                cur.active = true;
                cur.ter = pending_ter;
                pending_ter = false;
                std::memcpy(cur.name, rn, 3);
                cur.seq = seq;
                cur.chain = chain;
                cur.icode = icode;
            }

            const std::string_view atom = Field(p, len, 12, 4);
            int slot = -1;
            if (atom == "N") slot = 0;
            else if (atom == "CA") slot = 1;
            else if (atom == "C") slot = 2;

            if (slot >= 0) {
                const unsigned char rk = AltRank(alt);
                const unsigned char bit = static_cast<unsigned char>(1u << slot);
                if ((cur.mask & bit) == 0 || rk > cur.rank[slot]) {
                    float x = 0.0f;
                    float y = 0.0f;
                    float z = 0.0f;
                    if (ParseF(Field(p, len, 30, 8), x) &&
                        ParseF(Field(p, len, 38, 8), y) &&
                        ParseF(Field(p, len, 46, 8), z)) {
                        cur.xyz[slot][0] = x;
                        cur.xyz[slot][1] = y;
                        cur.xyz[slot][2] = z;
                        cur.mask = static_cast<unsigned char>(cur.mask | bit);
                        cur.rank[slot] = rk;
                    }
                }
            }
        }

        p = (nl != nullptr) ? nl + 1 : end;
    }

    commit(cur);
    return true;
}

void pbd_pool::Clear() {
    text_.clear();
    res_.clear();
}

std::size_t pbd_pool::capacity_bytes() const noexcept {
    return text_.capacity() + res_.capacity() * sizeof(Residue);
}

// ---------------------------------------------------------------- bin_pool

void bin_pool::Reset() {
    bytes_.clear();
    bytes_.resize(kHeaderBytes, 0u); // _len 占位, Save 时写回
    residues_ = 0;
}

void bin_pool::AppendF32(float v) {
    unsigned char raw[4];
    std::memcpy(raw, &v, sizeof raw);
    bytes_.insert(bytes_.end(), raw, raw + sizeof raw);
}

void bin_pool::Append(char kind, float cos_phi, float sin_phi, float cos_psi, float sin_psi) {
    bytes_.push_back(static_cast<unsigned char>(kind));
    AppendF32(cos_phi);
    AppendF32(sin_phi);
    AppendF32(cos_psi);
    AppendF32(sin_psi);
    ++residues_;
}

bool bin_pool::Save(const std::filesystem::path& out, std::string& err) const {
    err.clear();
    std::ofstream f(out, std::ios::binary | std::ios::trunc);
    if (!f) {
        err = "无法创建输出文件: " + Utf8FromPath(out);
        return false;
    }
    const std::int32_t len = static_cast<std::int32_t>(residues_);
    f.write(reinterpret_cast<const char*>(&len), sizeof len);
    if (bytes_.size() > kHeaderBytes) {
        f.write(reinterpret_cast<const char*>(bytes_.data() + kHeaderBytes),
                static_cast<std::streamsize>(bytes_.size() - kHeaderBytes));
    }
    f.flush();
    if (!f) {
        err = "写输出文件失败: " + Utf8FromPath(out);
        return false;
    }
    return true;
}

// ---------------------------------------------------------------- maker

void SimgleAnle_maker::Init(std::string bin_path) {
    bin_dir_ = PathFromUtf8(bin_path);
}

bool SimgleAnle_maker::Work(const std::filesystem::path& file,
                            std::uint32_t& residue_count,
                            std::uint64_t& out_bytes,
                            std::string& err) {
    residue_count = 0;
    out_bytes = 0;

    if (!pbd_.Load(file, err)) return false;
    if (!pbd_.Parse(err)) return false;

    const std::vector<Residue>& rs = pbd_.residues();
    bin_.Reset();
    for (std::size_t i = 0; i < rs.size(); ++i) {
        const Residue& r = rs[i];
        const bool full = (r.mask == 0x07u);

        float cos_phi = 0.0f;
        float sin_phi = 0.0f;
        float cos_psi = 0.0f;
        float sin_psi = 0.0f;

        if (full) {
            if (i > 0 && r.brk == 0 && (rs[i - 1].mask & 0x04u) != 0) {
                const float a = Dihedral(rs[i - 1].c, r.n, r.ca, r.c);
                cos_phi = std::cos(a);
                sin_phi = std::sin(a);
            }
            if (i + 1 < rs.size() && rs[i + 1].brk == 0 && (rs[i + 1].mask & 0x01u) != 0) {
                const float a = Dihedral(r.n, r.ca, r.c, rs[i + 1].n);
                cos_psi = std::cos(a);
                sin_psi = std::sin(a);
            }
        }

        bin_.Append(r.kind, cos_phi, sin_phi, cos_psi, sin_psi);
    }

    if (!bin_.Save(bin_dir_ / BinNameFor(file), err)) return false;

    residue_count = static_cast<std::uint32_t>(rs.size());
    out_bytes = static_cast<std::uint64_t>(kHeaderBytes + bin_.payload_bytes());
    return true;
}

std::size_t SimgleAnle_maker::capacity_bytes() const noexcept {
    return pbd_.capacity_bytes() + bin_.capacity_bytes();
}

} // namespace pbd
