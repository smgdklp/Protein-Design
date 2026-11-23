// von_dpmm / io.cpp — bin 读取、CSR、Bessel 表、二进制写出
#include "model.h"

#include <cstdio>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <filesystem>
#include <fstream>

namespace fs = std::filesystem;

namespace vd {

// ---------------------------------------------------------------- Bessel
// log( e^-k * I_m(k) ),直接级数 + logsumexp(精确,慢)
double logIe_exact(int m, double kappa) {
    if (kappa <= 1e-9) return (m == 0) ? 0.0 : -1e30;
    const double half = kappa * 0.5;
    const double lh = std::log(half);
    // j 上限:峰值在 j ~ kappa/2
    int jmax = int(kappa * 1.2) + 80 + 4 * m;
    if (jmax > 20000) jmax = 20000;
    double best = -1e300;
    std::vector<double> terms(jmax + 1);
    for (int j = 0; j <= jmax; ++j) {
        double t = (2.0 * j + m) * lh - std::lgamma(double(j) + 1.0)
                 - std::lgamma(double(j + m) + 1.0);
        terms[j] = t;
        if (t > best) best = t;
    }
    if (best < -1e299) return -1e30;
    double s = 0.0;
    for (int j = 0; j <= jmax; ++j) s += std::exp(terms[j] - best);
    return best + std::log(s) - kappa;
}

// kappa > 700 时 e^kappa 溢出 double,用 Debye 一致渐近(logIe 本身不溢出)
static double logIe_debye(int m, double kappa) {
    double mm = double(m);
    double r = std::sqrt(kappa * kappa + mm * mm);
    double asinh = std::log((mm + r) / kappa);      // asinh(m/k)
    return -0.5 * std::log(kTwoPi) - 0.5 * std::log(r) + (r - mm * asinh - kappa);
}

void build_bessel_table(BesselTable& t, int nk, double kmin, double kmax, int mmax) {
    t.nk = nk; t.mmax = mmax; t.kmin = kmin; t.kmax = kmax;
    t.dlogk = (std::log(kmax) - std::log(kmin)) / double(nk - 1);
    t.v.assign(size_t(nk) * size_t(mmax + 1), -1e30f);
    t.max_ref_err = 0.0;
    for (int i = 0; i < nk; ++i) {
        double k = kmin * std::exp(t.dlogk * double(i));
        for (int m = 0; m <= mmax; ++m) {
            double e = (k <= 700.0) ? logIe_exact(m, k) : logIe_debye(m, k);
            t.v[size_t(i) * size_t(mmax + 1) + m] = float(e);
            if (k <= 700.0 && m <= 8) {
                double r = logIe_exact(m, k);
                t.max_ref_err = std::max(t.max_ref_err, std::fabs(r - e));
            }
        }
    }
}

// vM sine 的 log 归一化常数(宿主参考版)
double logC_exact(double k1, double k2, double lam) {
    if (k1 <= 0) k1 = 1e-8;
    if (k2 <= 0) k2 = 1e-8;
    double t = lam * lam / (4.0 * k1 * k2);
    if (t > 1.0 / 4.0) t = 1.0 / 4.0;
    double lt = (t > 0) ? std::log(t) : -1e300;
    int M = int(8.0 * std::sqrt(std::max(k1, k2))) + 40;
    if (M > 220) M = 220;
    std::vector<double> terms(M + 1);
    double best = -1e300;
    for (int m = 0; m <= M; ++m) {
        double s = std::lgamma(2.0 * m + 1.0) - 2.0 * std::lgamma(double(m) + 1.0);
        if (m > 0) {
            if (lt < -1e299) { terms[m] = -1e300; continue; }
            s += m * lt;
        }
        s += logIe_exact(m, k1) + logIe_exact(m, k2);
        terms[m] = s;
        if (s > best) best = s;
    }
    double sum = 0.0;
    for (int m = 0; m <= M; ++m) sum += std::exp(terms[m] - best);
    double logS = best + std::log(sum);
    return -std::log(4.0 * kPi * kPi) - (k1 + k2) - logS;
}

// ---------------------------------------------------------------- 读 bin
namespace {

struct BinScanner {
    FILE* f = nullptr;
    std::vector<uint8_t> buf;
};

inline uint32_t rd_u32(const uint8_t* p) {
    uint32_t v; std::memcpy(&v, p, 4); return v;
}
inline float rd_f32(const uint8_t* p) {
    float v; std::memcpy(&v, p, 4); return v;
}

}  // namespace

bool load_points(const Config& cfg, PointSet& ps, std::string& err) {
    std::error_code ec;
    fs::path dir(cfg.bin_dir);
    if (!fs::exists(dir, ec)) { err = "bin dir missing: " + cfg.bin_dir; return false; }

    std::vector<fs::path> files;
    for (auto& e : fs::directory_iterator(dir, ec)) {
        if (!e.is_regular_file()) continue;
        if (e.path().extension() == ".bin") files.push_back(e.path());
    }
    if (files.empty()) { err = "no .bin under " + cfg.bin_dir; return false; }
    std::sort(files.begin(), files.end());   // 确定性

    const size_t expect = size_t(8'600'000);
    ps.ctx.reserve(expect);
    ps.cphi.reserve(expect); ps.sphi.reserve(expect);
    ps.cpsi.reserve(expect); ps.spsi.reserve(expect);
    ps.bx.reserve(expect);   ps.by.reserve(expect);

    std::vector<uint8_t> data;
    for (const fs::path& p : files) {
        std::ifstream in(p, std::ios::binary);
        if (!in) { err = "open failed: " + p.string(); return false; }
        in.seekg(0, std::ios::end);
        std::streamoff sz = in.tellg();
        in.seekg(0, std::ios::beg);
        data.resize(size_t(sz));
        if (sz > 0) in.read(reinterpret_cast<char*>(data.data()), sz);
        if (!in && sz > 0) { err = "read failed: " + p.string(); return false; }

        size_t off = 0;
        const size_t N = data.size();
        while (off + 4 <= N) {
            int32_t len = int32_t(rd_u32(&data[off]));
            off += 4;
            if (len <= 0) { err = fmt("bad len %d in %s", len, p.string().c_str()); return false; }
            size_t need = size_t(len) * 17;
            if (off + need > N) { err = fmt("truncated Que in %s", p.string().c_str()); return false; }

            ++ps.n_seg;
            ps.n_res += uint64_t(len);

            if (len >= 3) {
                // 为每个内部残基取三元组
                for (int32_t i = 1; i + 1 < len; ++i) {
                    const uint8_t* pc = &data[off + size_t(i - 1) * 17];
                    const uint8_t* cc = &data[off + size_t(i) * 17];
                    const uint8_t* nc = &data[off + size_t(i + 1) * 17];
                    int ia = aa_index(char(pc[0]));
                    int ib = aa_index(char(cc[0]));
                    int ic = aa_index(char(nc[0]));
                    if (ia < 0 || ib < 0 || ic < 0) { ++ps.drop_kine; continue; }

                    float cphi = rd_f32(cc + 1), sphi = rd_f32(cc + 5);
                    float cpsi = rd_f32(cc + 9), spsi = rd_f32(cc + 13);
                    if (!std::isfinite(cphi) || !std::isfinite(sphi) ||
                        !std::isfinite(cpsi) || !std::isfinite(spsi)) { ++ps.drop_nan; continue; }
                    double n1 = double(cphi) * cphi + double(sphi) * sphi;
                    double n2 = double(cpsi) * cpsi + double(spsi) * spsi;
                    if (std::fabs(n1 - 1.0) > 1e-6 || std::fabs(n2 - 1.0) > 1e-6) {
                        ++ps.drop_norm; continue;
                    }
                    ps.ctx.push_back(uint32_t(ctx_index(ia, ib, ic)));
                    ps.cphi.push_back(cphi); ps.sphi.push_back(sphi);
                    ps.cpsi.push_back(cpsi); ps.spsi.push_back(spsi);
                    ps.bx.push_back(uint8_t(bin_index(std::atan2(double(sphi), double(cphi)))));
                    ps.by.push_back(uint8_t(bin_index(std::atan2(double(spsi), double(cpsi)))));
                }
                ps.drop_edge += 2;
            } else {
                ps.drop_edge += uint64_t(len);
            }
            off += need;
        }
    }
    ps.n_points = ps.ctx.size();
    if (ps.n_points == 0) { err = "no usable points"; return false; }
    return true;
}

bool build_csr(const Config& cfg, PointSet& ps, std::string& err, CtxStat& st) {
    (void)cfg;
    ps.start.assign(size_t(kNctx) + 1, 0u);
    for (uint32_t c : ps.ctx) ps.start[size_t(c) + 1]++;
    for (int c = 0; c < kNctx; ++c) ps.start[size_t(c) + 1] += ps.start[size_t(c)];
    ps.order.assign(ps.n_points, 0u);

    std::vector<uint32_t> cursor(ps.start.begin(), ps.start.end() - 1);
    for (size_t i = 0; i < ps.n_points; ++i) {
        uint32_t c = ps.ctx[i];
        ps.order[cursor[c]++] = uint32_t(i);
    }
    st.n.assign(size_t(kNctx), 0u);
    for (int c = 0; c < kNctx; ++c) st.n[size_t(c)] = ps.start[size_t(c) + 1] - ps.start[size_t(c)];

    // 自检:单调、总数一致
    uint64_t sum = 0;
    for (int c = 0; c < kNctx; ++c) {
        if (ps.start[size_t(c) + 1] < ps.start[size_t(c)]) { err = "csr not monotone"; return false; }
        sum += st.n[size_t(c)];
    }
    if (sum != ps.n_points) { err = "csr total mismatch"; return false; }
    return true;
}

// ---------------------------------------------------------------- 写出
static bool write_box(const std::string& path, const char* magic,
                      const uint32_t* counts, uint64_t total, std::string& err) {
    FILE* f = fopen(path.c_str(), "wb");
    if (!f) { err = "open failed: " + path; return false; }
    BoxHeader h{};
    std::memcpy(h.magic, magic, 4);
    h.version = 1;
    h.bins = uint32_t(kBins);
    h.n_ctx = uint32_t(kNctx);
    h.n_cells = uint32_t(kCells);
    h.encoding = 1;
    h.reserved0 = 0;
    std::memcpy(h.aa, kAaOrder, 20);
    h.total = total;
    h.flags = 1;
    h.header_bytes = uint32_t(sizeof(BoxHeader));
    bool ok = fwrite(&h, sizeof(h), 1, f) == 1;

    size_t n = size_t(kNctx) * size_t(kCells);
    if (ok) ok = fwrite(counts, sizeof(uint32_t), n, f) == n;

    // 概率块:每个切片按质量归一
    std::vector<float> prob(n);
    for (int c = 0; c < kNctx && ok; ++c) {
        const uint32_t* src = counts + size_t(c) * kCells;
        double s = 0.0;
        for (int i = 0; i < kCells; ++i) s += double(src[i]);
        float inv = (s > 0) ? float(1.0 / s) : 0.0f;
        for (int i = 0; i < kCells; ++i) prob[size_t(c) * kCells + i] = float(src[i]) * inv;
    }
    if (ok) ok = fwrite(prob.data(), sizeof(float), n, f) == n;
    fclose(f);
    if (!ok) { err = "write failed: " + path; return false; }
    return true;
}

bool write_staic_box(const std::string& path, const uint32_t* counts, uint64_t total, std::string& err) {
    return write_box(path, "SBX1", counts, total, err);
}
bool write_von_box(const std::string& path, const uint32_t* counts, uint64_t total, std::string& err) {
    return write_box(path, "VBX1", counts, total, err);
}

bool write_params(const std::string& path, const ParamRec* recs, uint64_t total_points, std::string& err) {
    FILE* f = fopen(path.c_str(), "wb");
    if (!f) { err = "open failed: " + path; return false; }
    ParamHeader h{};
    std::memcpy(h.magic, "VPM1", 4);
    h.version = 2;
    h.bins = uint32_t(kBins);
    h.n_ctx = uint32_t(kNctx);
    h.k_max = uint32_t(kKMax);
    h.n_aa = uint32_t(kNaa);
    std::memcpy(h.aa, kAaOrder, 20);
    h.total_points = total_points;
    h.flags = 1;
    h.rec_bytes = uint32_t(sizeof(ParamRec));
    bool ok = fwrite(&h, sizeof(h), 1, f) == 1;
    if (ok) ok = fwrite(recs, sizeof(ParamRec), size_t(kNctx), f) == size_t(kNctx);
    fclose(f);
    if (!ok) { err = "write failed: " + path; return false; }
    return true;
}

bool write_pred_grid(const std::string& path, const float* grid, std::string& err) {
    FILE* f = fopen(path.c_str(), "wb");
    if (!f) { err = "open failed: " + path; return false; }
    struct H { char magic[4]; uint32_t version, bins, n_ctx, n_cells, enc; uint64_t total; uint32_t hb, rsv; };
    H h{};
    std::memcpy(h.magic, "VPG1", 4);
    h.version = 1; h.bins = uint32_t(kBins); h.n_ctx = uint32_t(kNctx);
    h.n_cells = uint32_t(kCells); h.enc = 0; h.total = 0;
    h.hb = uint32_t(sizeof(H)); h.rsv = 0;
    bool ok = fwrite(&h, sizeof(h), 1, f) == 1;
    size_t n = size_t(kNctx) * size_t(kCells);
    if (ok) ok = fwrite(grid, sizeof(float), n, f) == n;
    fclose(f);
    if (!ok) { err = "write failed: " + path; return false; }
    return true;
}

}  // namespace vd
