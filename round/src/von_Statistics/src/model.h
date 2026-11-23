// von_dpmm / model.h — 常量、索引、配置、二进制格式、宿主侧声明
#pragma once

#include <cstdint>
#include <cstddef>
#include <cmath>
#include <string>
#include <vector>
#include <array>

namespace vd {

// ---------------------------------------------------------------- 基本常量
constexpr int  kBins    = 36;              // 每轴 36 格 = 10°
constexpr int  kCells   = kBins * kBins;   // 1296
constexpr int  kNaa     = 20;
constexpr int  kNctx    = kNaa * kNaa * kNaa;   // 8000
constexpr int  kKMax    = 8;               // 压缩后的最大成分数
constexpr int  kKCap    = 24;              // MCMC 内最大占用簇数

constexpr double kPi    = 3.14159265358979323846;
constexpr double kTwoPi = 6.28318530717958647692;

inline const char* kAaOrder = "ACDEFGHIKLMNPQRSTVWY";

inline int aa_index(char c) {
    for (int i = 0; i < kNaa; ++i) if (kAaOrder[i] == c) return i;
    return -1;
}
inline int ctx_index(int a, int b, int c) { return (a * kNaa + b) * kNaa + c; }

// 角度 -> 格号:[-pi,pi) 均分 36 格,pi 归到最后一格
inline int bin_index(double ang) {
    double t = (ang + kPi) / kTwoPi * double(kBins);
    int k = int(std::floor(t));
    if (k < 0) k = 0;
    if (k >= kBins) k = kBins - 1;
    return k;
}
inline double bin_center(int k) { return -kPi + (double(k) + 0.5) * (kTwoPi / double(kBins)); }
constexpr double kCellArea = (kTwoPi / double(kBins)) * (kTwoPi / double(kBins));

// ---------------------------------------------------------------- 配置
struct Config {
    // 数据源
    std::string bin_dir   = "H:\\DATABACE\\init";
    std::string manifest  = "H:\\DATABACE\\init\\manifest.tsv";
    std::string out_root  = "H:\\DATABACE\\thr_result";

    // 分箱
    int bins = kBins;

    // MCMC 种子
    uint64_t seed = 20240601ull;

    // 分级阈值(样本数)
    int tier_a_n = 2000;
    int tier_b_n = 200;

    // 每级迭代预算 {iters, burn, thin, chains}
    int a_iters = 6000, a_burn = 1000, a_thin = 10, a_chains = 2;
    int b_iters = 4000, b_burn = 1000, b_thin = 5,  b_chains = 1;
    int c_iters = 500,  c_burn = 100,  c_thin = 5,  c_chains = 1;

    // DP 质量参数 / 先验
    float tau0   = 1.0f;
    float k10    = 0.10132118f;   // 1/pi^2
    float k20    = 0.10132118f;
    float lam0   = 0.0f;
    int   wish_nu0 = 2;
    float wish_b0  = 0.25f;       // 0.5^2 * I

    // 采样与压缩
    int      k_aux      = 3;      // Algorithm 8 辅助簇数
    int      k_max_comp = 8;      // 压缩后最大成分
    uint64_t m_samples  = 200000; // von_box 每上下文采样数

    // 运行控制
    bool cpu_only   = false;
    int  max_ctx    = 0;          // >0 时只跑前 N 个上下文(调试)
    bool skip_fit   = false;
    bool skip_sample= false;
};

// ---------------------------------------------------------------- 数据容器(宿主)
struct PointSet {
    std::vector<uint32_t> ctx;   // 上下文 id
    std::vector<uint32_t> order; // 按 ctx 排序后的点索引
    std::vector<uint32_t> start; // kNctx+1
    std::vector<float>    cphi, sphi, cpsi, spsi;
    std::vector<uint8_t>  bx, by;

    uint64_t n_points  = 0;
    uint64_t n_res     = 0;   // 扫描到的残基数
    uint64_t n_seg     = 0;   // 段数
    uint64_t drop_edge = 0;   // 边界(首尾残基)
    uint64_t drop_nan  = 0;   // 角度非法
    uint64_t drop_norm = 0;   // cos^2+sin^2 超差
    uint64_t drop_kine = 0;   // 非 20 标准

    size_t size() const { return size_t(n_points); }
};

struct CtxStat {
    std::vector<uint32_t> n;   // 每个 ctx 的点数
};

// ---------------------------------------------------------------- Bessel 表
// 存放 log(e^-k * I_m(k)),用于 vM sine 归一化常数 C 的快速求值。
struct BesselTable {
    int    nk   = 0;      // kappa 网格点数
    int    mmax = 0;      // 最大阶 m
    double kmin = 1e-3;
    double kmax = 6000.0;
    double dlogk = 0.0;
    std::vector<float> v; // [nk][mmax+1]
    double max_ref_err = 0.0;   // 建表时与直算 logIe 的最大偏差(自检)
};

// 单点精确值(直接级数,慢,只用于自检与建表)
double logIe_exact(int m, double kappa);
void   build_bessel_table(BesselTable& t, int nk = 2048, double kmin = 1e-3,
                          double kmax = 6000.0, int mmax = 220);

// vM sine 对数归一化(宿主参考实现,用于自检)
double logC_exact(double k1, double k2, double lam);

// ---------------------------------------------------------------- 二进制格式
struct BoxHeader {
    char     magic[4];      // "SBX1" / "VBX1"
    uint32_t version;
    uint32_t bins;
    uint32_t n_ctx;
    uint32_t n_cells;
    uint32_t encoding;      // 1 = counts(u32) + probs(f32)
    uint32_t reserved0;
    char     aa[20];
    uint64_t total;         // 总样本/点数
    uint32_t flags;         // bit0: 每片已归一
    uint32_t header_bytes;  // 64
};
static_assert(sizeof(BoxHeader) == 64, "BoxHeader must be 64 bytes");

struct ParamRec {
    uint32_t tier;          // 0=D(空),1=C,2=B,3=A
    uint32_t n_obs;
    uint32_t k_eff;
    uint32_t n_samples;
    float    k_post_mean;
    float    k_post_sd;
    float    loglik_mean;
    float    reserved_f;
    // [k][8]: w, mu, nu, kappa1, kappa2, lambda, logC, reserved
    float    comp[kKMax * 8];
};
static_assert(sizeof(ParamRec) == 288, "ParamRec must be 288 bytes");

struct ParamHeader {
    char     magic[4];      // "VPM1"
    uint32_t version;
    uint32_t bins;
    uint32_t n_ctx;
    uint32_t k_max;
    uint32_t n_aa;
    char     aa[20];
    uint64_t total_points;
    uint32_t flags;
    uint32_t rec_bytes;
    uint64_t reserved;
};
static_assert(sizeof(ParamHeader) == 72, "ParamHeader must be 72 bytes");

// ---------------------------------------------------------------- 接口
bool load_config(const std::string& path, Config& cfg, std::string& err);

bool load_points(const Config& cfg, PointSet& ps, std::string& err);
bool build_csr(const Config& cfg, PointSet& ps, std::string& err, CtxStat& st);

bool write_staic_box(const std::string& path, const uint32_t* counts, uint64_t total, std::string& err);
bool write_von_box  (const std::string& path, const uint32_t* counts, uint64_t total, std::string& err);
bool write_params   (const std::string& path, const ParamRec* recs, uint64_t total_points, std::string& err);
bool write_pred_grid(const std::string& path, const float* grid, std::string& err);

std::string fmt(const char* f, ...);

}  // namespace vd
