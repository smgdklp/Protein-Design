// von_dpmm / cpu.cpp — 与 GPU 核同算法的宿主实现(参考实现 + CUDA 不可用时的退路)
#include "model.h"
#include "gpu.h"

#include <cstdio>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <thread>
#include <vector>

namespace vd {
namespace {

struct Rng { uint64_t s; };

inline uint64_t mix64(uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}
inline Rng make_rng(uint64_t seed, uint32_t ctx, uint32_t chain, uint32_t iter, uint32_t salt) {
    Rng r;
    uint64_t s = seed
               ^ (uint64_t(ctx + 1) * 0x9E3779B97F4A7C15ULL)
               ^ (uint64_t(chain + 1) * 0xD1B54A32D192ED03ULL)
               ^ (uint64_t(iter + 1) * 0x94D049BB133111EBULL)
               ^ (uint64_t(salt + 1) * 0xA24BAED4963EE407ULL);
    r.s = mix64(s);
    return r;
}
inline double u01(Rng& r) {
    r.s = mix64(r.s);
    return double((r.s >> 11) & 0x1FFFFFFFFFFFFFULL) * (1.0 / 9007199254740992.0);
}
inline double nrm(Rng& r) {
    double u1 = std::max(u01(r), 1e-12), u2 = u01(r);
    return std::sqrt(-2.0 * std::log(u1)) * std::cos(kTwoPi * u2);
}

// 表插值 log(e^-k I_m(k))
inline double logIe_tab(const BesselTable& tab, int m, double k) {
    if (m > tab.mmax) return -1e30;
    double lk = std::log(std::max(k, 1e-9));
    double t = (lk - std::log(tab.kmin)) / tab.dlogk;
    if (t < 0) t = 0;
    if (t > double(tab.nk - 1)) t = double(tab.nk - 1);
    int i0 = int(t);
    int i1 = std::min(i0 + 1, tab.nk - 1);
    double f = t - double(i0);
    const float* r0 = tab.v.data() + size_t(i0) * size_t(tab.mmax + 1) + m;
    const float* r1 = tab.v.data() + size_t(i1) * size_t(tab.mmax + 1) + m;
    return double(r0[0]) * (1.0 - f) + double(r1[0]) * f;
}

// log C(2m,m) 只与 m 有关:预先算好。否则每个 MH 提议都要花 2*M 次 lgamma,
// 单次 logC 能吃掉几十微秒,是 CPU 后端最大的隐性瓶颈。
struct LogBinom {
    double v[256];
    LogBinom() {
        for (int m = 0; m < 256; ++m)
            v[m] = std::lgamma(2.0 * m + 1.0) - 2.0 * std::lgamma(double(m) + 1.0);
    }
};
inline const double* logbinom_table() {
    static const LogBinom t;
    return t.v;
}

inline double logC_tab(const BesselTable& tab, double k1, double k2, double lam) {
    k1 = std::max(k1, 1e-7);
    k2 = std::max(k2, 1e-7);
    double t = lam * lam / (4.0 * k1 * k2);
    if (t > 0.2499999) t = 0.2499999;
    double lt = (t > 1e-30) ? std::log(t) : -1e300;
    int M = int(8.0 * std::sqrt(std::max(k1, k2))) + 40;
    if (M > tab.mmax) M = tab.mmax;
    if (M > 255) M = 255;
    const double* lb = logbinom_table();
    double best = -1e300, sum = 0.0;
    for (int m = 0; m <= M; ++m) {
        double s = lb[m];
        if (m > 0) {
            if (lt < -1e299) break;
            s += double(m) * lt;
        }
        s += logIe_tab(tab, m, k1) + logIe_tab(tab, m, k2);
        if (s > best) { sum = sum * std::exp(best - s) + 1.0; best = s; }
        else sum += std::exp(s - best);
    }
    return -std::log(4.0 * kPi * kPi) - (k1 + k2) - (best + std::log(std::max(sum, 1e-300)));
}

// 单变量 von Mises 采样:查表反 CDF。
// (Best-Fisher 的拒绝采样版本我第一版写错了,在 kappa>=10 时 100% 退化到 mu±pi;
//  这里换成数值反 CDF:kappa 取对数网格 256 档,半圆 [0,pi] 上 2048 个节点,
//  对称性决定了只需半圆 + 随机符号。分辨率 0.09°,远小于 10° 的格子。)
struct VmTable {
    static const int NK = 256;
    static const int NT = 2048;
    double lkmin, lkmax;
    std::vector<std::vector<float>> cdf;
    VmTable() {
        lkmin = std::log(1e-3);
        lkmax = std::log(2000.0);
        cdf.resize(NK);
        for (int i = 0; i < NK; ++i) {
            double k = std::exp(lkmin + (lkmax - lkmin) * double(i) / double(NK - 1));
            cdf[i].resize(NT + 1);
            double sum = 0.0;
            cdf[i][0] = 0.f;
            for (int j = 1; j <= NT; ++j) {
                double t = (double(j) - 0.5) * kPi / double(NT);
                sum += std::exp(k * std::cos(t) - k);      // 减去 k 防止溢出
                cdf[i][j] = float(sum);
            }
            if (sum > 0) for (int j = 0; j <= NT; ++j) cdf[i][j] = float(double(cdf[i][j]) / sum);
            cdf[i][NT] = 1.f;
        }
    }
    inline double sample_abs(double kappa, Rng& r) const {
        double lk = std::log(kappa);
        double fi = (lk - lkmin) / (lkmax - lkmin) * double(NK - 1);
        if (fi < 0) fi = 0;
        if (fi > double(NK - 1)) fi = double(NK - 1);
        int i = int(fi);
        const std::vector<float>& c = cdf[size_t(i)];
        double u = u01(r);
        int lo = 0, hi = NT;
        while (lo + 1 < hi) { int mid = (lo + hi) >> 1; if (double(c[size_t(mid)]) < u) lo = mid; else hi = mid; }
        double c0 = c[size_t(lo)], c1 = c[size_t(hi)];
        double frac = (c1 > c0) ? (u - c0) / (c1 - c0) : 0.0;
        return (double(lo) + frac) * kPi / double(NT);
    }
};

inline double sample_vm(double mu, double kappa, Rng& r) {
    kappa = std::fabs(kappa);
    if (!(kappa > 1e-6)) return mu + kTwoPi * u01(r);
    static const VmTable tab;
    double d = (kappa > 2000.0) ? tab.sample_abs(2000.0, r) : tab.sample_abs(kappa, r);
    return mu + ((u01(r) < 0.5) ? -d : d);
}

inline double invA(double R) {
    R = std::min(std::max(R, 0.0), 0.999999);
    if (R < 0.53) return 2.0 * R + R * R * R + 5.0 * R * R * R * R * R / 6.0;
    if (R < 0.85) return -0.4 + 1.39 * R + 0.43 / (1.0 - R);
    return 1.0 / std::max(R * R * R - 4.0 * R * R + 3.0 * R, 1e-6);
}

struct FitIn {
    const PointSet* ps;
    const CtxStat*  st;
    const Config*   cfg;
    const BesselTable* tab;
    std::vector<ParamRec>* recs;
    std::vector<float>*    pred;
};

struct CtxItem { int ctx; int tier; int n; };

// ---------------------------------------------------------------- 单上下文 MCMC
void run_one_ctx(const FitIn& F, int ctx, int tier, std::vector<float>& gridOut, ParamRec& rec, std::vector<std::array<double, 6>>& seeds) {
    const Config& cfg = *F.cfg;
    const BesselTable& tab = *F.tab;
    const PointSet& ps = *F.ps;
    const uint32_t p0 = ps.start[size_t(ctx)];
    const int n = int(ps.start[size_t(ctx) + 1] - p0);

    int iters, burn, thin, chains;
    if (tier >= 3)      { iters = cfg.a_iters; burn = cfg.a_burn; thin = cfg.a_thin; chains = cfg.a_chains; }
    else if (tier == 2) { iters = cfg.b_iters; burn = cfg.b_burn; thin = cfg.b_thin; chains = cfg.b_chains; }
    else                { iters = cfg.c_iters; burn = cfg.c_burn; thin = cfg.c_thin; chains = cfg.c_chains; }

    const int K10 = 0, K20 = 0;   // 未用
    (void)K10; (void)K20;
    const double c10 = cfg.k10, c20 = cfg.k20, l0 = cfg.lam0, tau0 = cfg.tau0;
    const int KAUX = std::min(cfg.k_aux, 8);

    std::vector<int> lab(size_t(n), 0);
    double s_mu[kKCap], s_nu[kKCap], s_k1[kKCap], s_k2[kKCap], s_lam[kKCap], s_lc[kKCap];
    int    s_cnt[kKCap];
    double a_mu[8], a_nu[8], a_k1[8], a_k2[8], a_lam[8], a_lc[8];

    std::vector<float> grid(size_t(kCells), 0.f);
    int K_end = 1;
    double km_sum = 0.0, km2_sum = 0.0;
    int    nkeep = 0;

    for (int chain = 0; chain < chains; ++chain) {
        int K = 1;
        s_cnt[0] = n;
        s_mu[0] = 0; s_nu[0] = 0;
        s_k1[0] = 1; s_k2[0] = 1; s_lam[0] = 0;
        s_lc[0] = logC_tab(tab, 1, 1, 0);
        std::fill(lab.begin(), lab.end(), 0);

        double acc[3] = {0.0, 0.0, 0.0};   // 未用

        for (int it = 0; it < iters; ++it) {
            // ---- 辅助簇
            for (int j = 0; j < KAUX; ++j) {
                Rng r = make_rng(cfg.seed, uint32_t(ctx), uint32_t(chain), uint32_t(it), uint32_t(j) + 17u);
                double mu = sample_vm(0.0, c10, r);
                double nu = sample_vm(0.0, c20, r);
                double o11 = cfg.wish_b0, o22 = cfg.wish_b0, o12 = 0.0;
                {
                    int nu_ = std::max(1, cfg.wish_nu0);
                    double c2 = 0.0;
                    for (int q = 0; q < nu_; ++q) { double z = nrm(r); c2 += z * z; }
                    c2 = std::max(c2, 1e-9);
                    double l11 = std::sqrt(std::max(cfg.wish_b0 * c2, 1e-12));
                    double z2 = nrm(r);
                    double l21 = 0.0 * std::sqrt(c2) / l11 + std::sqrt(std::max(double(cfg.wish_b0), 1e-9)) * z2;
                    double cond = std::max(double(cfg.wish_b0), 1e-9);
                    double l22 = std::sqrt(std::max(cond * std::max(c2 - z2 * z2, 1e-6), 1e-6));
                    o11 = l11 * l11; o12 = l11 * l21; o22 = l21 * l21 + l22 * l22;
                }
                double k1 = o11, k2 = std::max(o22, 1e-4), lam = -o12;
                if (lam * lam > 0.98 * k1 * k2) lam = std::sqrt(0.98 * k1 * k2) * (lam >= 0 ? 1.0 : -1.0);
                a_mu[j] = mu; a_nu[j] = nu; a_k1[j] = k1; a_k2[j] = k2; a_lam[j] = lam;
                a_lc[j] = logC_tab(tab, k1, k2, lam);
            }

            const double lden = std::log(double(n) + tau0);
            const double lw_aux = std::log(tau0 / double(KAUX)) - lden;
            double lw_real[kKCap];
            for (int k = 0; k < K; ++k) lw_real[k] = std::log(double(s_cnt[k])) - lden;

            // 预处理每个成分的 cos/sin(点内循环不再做三角函数)
            double cmu[32], smu[32], cnu[32], snu[32];
            for (int k = 0; k < K; ++k) {
                cmu[k] = std::cos(s_mu[k]); smu[k] = std::sin(s_mu[k]);
                cnu[k] = std::cos(s_nu[k]); snu[k] = std::sin(s_nu[k]);
            }
            for (int j = 0; j < KAUX; ++j) {
                int k = K + j;
                cmu[k] = std::cos(a_mu[j]); smu[k] = std::sin(a_mu[j]);
                cnu[k] = std::cos(a_nu[j]); snu[k] = std::sin(a_nu[j]);
            }

            // ---- 分配
            for (int t = 0; t < n; ++t) {
                uint32_t gp = p0 + uint32_t(t);
                double cph = ps.cphi[gp], sph = ps.sphi[gp], cps = ps.cpsi[gp], sps = ps.spsi[gp];
                double sc[32];
                double best = -3e38;
                int ncomp = K + KAUX;
                for (int k = 0; k < ncomp; ++k) {
                    double k1, k2, lam, lc, lw;
                    if (k < K) { k1 = s_k1[k]; k2 = s_k2[k]; lam = s_lam[k]; lc = s_lc[k]; lw = lw_real[k]; }
                    else { int j = k - K; k1 = a_k1[j]; k2 = a_k2[j]; lam = a_lam[j]; lc = a_lc[j]; lw = lw_aux; }
                    double c1 = cph * cmu[k] + sph * smu[k], q1 = sph * cmu[k] - cph * smu[k];
                    double c2 = cps * cnu[k] + sps * snu[k], q2 = sps * cnu[k] - cps * snu[k];
                    double ld = lc + k1 * c1 + k2 * c2 + lam * q1 * q2 + lw;
                    sc[k] = ld;
                    if (ld > best) best = ld;
                }
                double sum = 0.0;
                for (int k = 0; k < ncomp; ++k) { sc[k] = double(expf(float(sc[k] - best))); sum += sc[k]; }
                Rng r = make_rng(cfg.seed, uint32_t(ctx), uint32_t(chain), uint32_t(it), uint32_t(t) * 3u + 1013u);
                double u = u01(r) * sum, a = 0.0;
                int pick = ncomp - 1;
                for (int k = 0; k < ncomp; ++k) { a += sc[k]; if (u <= a) { pick = k; break; } }
                lab[size_t(t)] = pick;
            }

            // ---- 落簇(辅助簇 -> 真实簇)
            {
                const int Kold = K;
                for (int j = 0; j < KAUX; ++j) {
                    int c = 0;
                    for (int t = 0; t < n; ++t) if (lab[size_t(t)] == Kold + j) ++c;
                    if (c > 0) {
                        if (K < kKCap) {
                            int k = K++;
                            s_mu[k] = a_mu[j]; s_nu[k] = a_nu[j];
                            s_k1[k] = a_k1[j]; s_k2[k] = a_k2[j]; s_lam[k] = a_lam[j]; s_lc[k] = a_lc[j];
                            s_cnt[k] = c;
                            for (int t = 0; t < n; ++t) if (lab[size_t(t)] == Kold + j) lab[size_t(t)] = k;
                        } else {
                            for (int t = 0; t < n; ++t) if (lab[size_t(t)] == Kold + j) lab[size_t(t)] = 0;
                        }
                    }
                }
            }
            for (int k = 0; k < K; ++k) s_cnt[k] = 0;
            for (int t = 0; t < n; ++t) { int k = lab[size_t(t)]; if (k >= 0 && k < K) s_cnt[k]++; }

            // ---- 回收空簇(K 只增不减会一路顶到 kKCap,导致过分裂 + 变慢)
            {
                int mp[kKCap];
                int nk2 = 0;
                for (int k = 0; k < K; ++k) mp[k] = (s_cnt[k] > 0) ? nk2++ : -1;
                if (nk2 < K) {
                    for (int k = 0; k < K; ++k) {
                        int d = mp[k];
                        if (d < 0 || d == k) continue;
                        s_mu[d] = s_mu[k]; s_nu[d] = s_nu[k];
                        s_k1[d] = s_k1[k]; s_k2[d] = s_k2[k];
                        s_lam[d] = s_lam[k]; s_lc[d] = s_lc[k];
                        s_cnt[d] = s_cnt[k];
                    }
                    for (int t = 0; t < n; ++t) {
                        int k = lab[size_t(t)];
                        int d = (k >= 0 && k < K) ? mp[k] : -1;
                        lab[size_t(t)] = (d >= 0) ? d : 0;
                    }
                    if (nk2 == 0) { s_cnt[0] = n; s_mu[0] = 0; s_nu[0] = 0; s_k1[0] = 1; s_k2[0] = 1; s_lam[0] = 0; s_lc[0] = logC_tab(tab, 1, 1, 0); }
                    K = nk2 > 0 ? nk2 : 1;
                }
            }

            // ---- μ|ν
            double sA[kKCap] = {0}, sB[kKCap] = {0};
            double csin[kKCap], ccos[kKCap];
            for (int k = 0; k < K; ++k) { csin[k] = std::sin(s_nu[k]); ccos[k] = std::cos(s_nu[k]); }
            for (int t = 0; t < n; ++t) {
                int k = lab[size_t(t)];
                if (k < 0 || k >= K) continue;
                uint32_t gp = p0 + uint32_t(t);
                double cph = ps.cphi[gp], sph = ps.sphi[gp], cps = ps.cpsi[gp], sps = ps.spsi[gp];
                double sn = csin[k], cn = ccos[k];
                double s2 = sps * cn - cps * sn;
                sA[k] += s_k1[k] * cph + s_lam[k] * s2 * sph;
                sB[k] += s_k1[k] * sph - s_lam[k] * s2 * cph;
            }
            for (int k = 0; k < K; ++k) {
                double Aa = sA[k] + c10;
                double Bb = sB[k] - l0 * std::sin(s_nu[k]);
                double kap = std::sqrt(Aa * Aa + Bb * Bb);
                Rng r = make_rng(cfg.seed, uint32_t(ctx), uint32_t(chain), uint32_t(it), uint32_t(k) + 7001u);
                s_mu[k] = sample_vm(std::atan2(Bb, Aa), kap, r);
            }
            // ---- ν|μ
            double sC[kKCap] = {0}, sD[kKCap] = {0};
            for (int k = 0; k < K; ++k) { csin[k] = std::sin(s_mu[k]); ccos[k] = std::cos(s_mu[k]); }
            for (int t = 0; t < n; ++t) {
                int k = lab[size_t(t)];
                if (k < 0 || k >= K) continue;
                uint32_t gp = p0 + uint32_t(t);
                double cph = ps.cphi[gp], sph = ps.sphi[gp], cps = ps.cpsi[gp], sps = ps.spsi[gp];
                double sm = csin[k], cm = ccos[k];
                double s1 = sph * cm - cph * sm;
                sC[k] += s_k2[k] * cps + s_lam[k] * s1 * sps;
                sD[k] += s_k2[k] * sps - s_lam[k] * s1 * cps;
            }
            for (int k = 0; k < K; ++k) {
                double Aa = sC[k] + c20;
                double Bb = sD[k] - l0 * std::sin(s_mu[k]);
                double kap = std::sqrt(Aa * Aa + Bb * Bb);
                Rng r = make_rng(cfg.seed, uint32_t(ctx), uint32_t(chain), uint32_t(it), uint32_t(k) + 9001u);
                s_nu[k] = sample_vm(std::atan2(Bb, Aa), kap, r);
            }
            // ---- 精度 MH
            {
                double S1[kKCap] = {0}, S2[kKCap] = {0}, S3[kKCap] = {0};
                for (int t = 0; t < n; ++t) {
                    int k = lab[size_t(t)];
                    if (k < 0 || k >= K) continue;
                    uint32_t gp = p0 + uint32_t(t);
                    double cph = ps.cphi[gp], sph = ps.sphi[gp], cps = ps.cpsi[gp], sps = ps.spsi[gp];
                    double cm = std::cos(s_mu[k]), sm = std::sin(s_mu[k]);
                    double cn = std::cos(s_nu[k]), sn = std::sin(s_nu[k]);
                    double c1 = cph * cm + sph * sm, q1 = sph * cm - cph * sm;
                    double c2 = cps * cn + sps * sn, q2 = sps * cn - cps * sn;
                    S1[k] += c1; S2[k] += c2; S3[k] += q1 * q2;
                }
                for (int k = 0; k < K; ++k) {
                    double k1 = s_k1[k], k2 = s_k2[k], lam = s_lam[k], lc0 = s_lc[k];
                    double cnt = double(s_cnt[k]);
                    double lp0 = k1 * S1[k] + k2 * S2[k] + lam * S3[k] + cnt * lc0
                               - 0.5 * std::log(std::max(k1 * k2 - lam * lam, 1e-12))
                               - (0.5 / cfg.wish_b0) * (k1 + k2);
                    // 每轮做 3 次随机游走 MH:精度参数必须能快速追上数据,
                    // 否则新簇停在"近均匀"处,混合永远尖不起来(过平滑的根因)。
                    Rng r = make_rng(cfg.seed, uint32_t(ctx), uint32_t(chain), uint32_t(it), uint32_t(k) + 12007u);
                    for (int rep = 0; rep < 3; ++rep) {
                        double rho = lam / std::sqrt(std::max(k1 * k2, 1e-12));
                        double z = std::atanh(std::min(std::max(rho, -0.999), 0.999));
                        double sc_ = 0.55;
                        double nk1 = k1 * std::exp(sc_ * nrm(r));
                        double nk2 = k2 * std::exp(sc_ * nrm(r));
                        double nz = z + sc_ * nrm(r);
                        double nrho = std::tanh(std::min(std::max(nz, -3.0), 3.0));
                        double nlam = nrho * std::sqrt(nk1 * nk2);
                        double nlc = logC_tab(tab, nk1, nk2, nlam);
                        double lp1 = nk1 * S1[k] + nk2 * S2[k] + nlam * S3[k] + cnt * nlc
                                   - 0.5 * std::log(std::max(nk1 * nk2 - nlam * nlam, 1e-12))
                                   - (0.5 / cfg.wish_b0) * (nk1 + nk2);
                        double u = u01(r);
                        if (std::log(std::max(u, 1e-12)) < (lp1 - lp0)) {
                            k1 = nk1; k2 = nk2; lam = nlam; lc0 = nlc;
                            lp0 = lp1;
                            s_k1[k] = nk1; s_k2[k] = nk2; s_lam[k] = nlam; s_lc[k] = nlc;
                        }
                    }
                }
            }

            // ---- 网格累积
            if (it >= burn && ((it - burn) % thin) == 0) {
                double wsum = 0.0;
                for (int k = 0; k < K; ++k) wsum += double(s_cnt[k]);
                double wnew = tau0 / (wsum + tau0);
                double lnew = logC_tab(tab, c10, c20, l0);
                for (int c = 0; c < kCells; ++c) {
                    int ix = c / kBins, iy = c % kBins;
                    double ph = bin_center(ix), qs = bin_center(iy);
                    double cph = std::cos(ph), sph = std::sin(ph), cps = std::cos(qs), sps = std::sin(qs);
                    double v = 0.0;
                    for (int k = 0; k < K; ++k) {
                        double cm = std::cos(s_mu[k]), sm = std::sin(s_mu[k]);
                        double cn = std::cos(s_nu[k]), sn = std::sin(s_nu[k]);
                        double c1 = cph * cm + sph * sm, q1 = sph * cm - cph * sm;
                        double c2 = cps * cn + sps * sn, q2 = sps * cn - cps * sn;
                        double w = double(s_cnt[k]) / (wsum + tau0);
                        v += w * std::exp(s_lc[k] + s_k1[k] * c1 + s_k2[k] * c2 + s_lam[k] * q1 * q2);
                    }
                    {
                        double c1 = cph, q1 = sph;
                        double c2 = cps, q2 = sps;
                        v += wnew * std::exp(lnew + c10 * c1 + c20 * c2 + l0 * q1 * q2);
                    }
                    grid[size_t(c)] += float(v);
                }
                km_sum += double(K); km2_sum += double(K) * double(K); ++nkeep;
            }
        }
        K_end = K;
    }

    gridOut = grid;
    seeds.clear();
    {
        std::vector<int> ord;
        ord.resize(size_t(K_end));
        for (int k = 0; k < K_end; ++k) ord[size_t(k)] = k;
        std::sort(ord.begin(), ord.end(), [&](int a, int b) { return s_cnt[a] > s_cnt[b]; });
        for (int i = 0; i < K_end && i < kKMax; ++i) {
            int k = ord[size_t(i)];
            seeds.push_back({double(s_cnt[k]) / (double(n) + tau0), s_mu[k], s_nu[k],
                             s_k1[k], s_k2[k], s_lam[k]});
        }
    }
    if (const char* dbg = std::getenv("VD_DEBUG_CTX")) {
        if (std::atoi(dbg) == ctx) {
            std::printf("\n[debug] ctx=%d n=%d K=%d nkeep=%d kpm=%.2f\n", ctx, n, K_end, nkeep,
                        nkeep > 0 ? km_sum / nkeep : 0.0);
            for (int k = 0; k < K_end; ++k) {
                std::printf("   k=%d cnt=%d mu=%+.3f nu=%+.3f k1=%8.2f k2=%8.2f lam=%+8.2f lc=%+.3f\n",
                            k, s_cnt[k], s_mu[k], s_nu[k], s_k1[k], s_k2[k], s_lam[k], s_lc[k]);
            }
            std::fflush(stdout);
        }
    }
    rec.tier = uint32_t(tier);
    rec.n_obs = uint32_t(n);
    rec.n_samples = uint32_t(nkeep);
    rec.k_eff = 0;
    double m = nkeep > 0 ? km_sum / nkeep : 0.0;
    double m2 = nkeep > 0 ? km2_sum / nkeep : 0.0;
    rec.k_post_mean = float(m);
    rec.k_post_sd = float(std::sqrt(std::max(m2 - m * m, 0.0)));
    rec.loglik_mean = 0.f;
    rec.reserved_f = 0.f;
    for (int k = 0; k < kKMax; ++k) {
        float* d = rec.comp + k * 8;
        d[0] = 0.f; d[1] = 0.f; d[2] = 0.f; d[3] = 1.f; d[4] = 1.f; d[5] = 0.f; d[6] = 0.f; d[7] = 0.f;
    }
}

// ---------------------------------------------------------------- 压缩
// seeds:MCMC 末态的簇参数(按权重降序),每簇 6 个数 {w,mu,nu,k1,k2,lam}
void compress_ctx(const BesselTable& tab, int n_obs, const std::vector<float>& grid,
                  int kmax, int em_iters, ParamRec& out,
                  const std::vector<std::array<double, 6>>& seeds) {
    double mass[kCells];
    double tot = 0.0;
    for (int c = 0; c < kCells; ++c) tot += grid[size_t(c)];
    if (tot <= 0) tot = 1.0;
    for (int c = 0; c < kCells; ++c) mass[c] = double(grid[size_t(c)]) / tot;

    double cph[kBins], sph[kBins];
    for (int k = 0; k < kBins; ++k) { cph[k] = std::cos(bin_center(k)); sph[k] = std::sin(bin_center(k)); }

    double best_bic = 3e38;
    ParamRec best{};
    best.tier = out.tier; best.n_obs = out.n_obs; best.n_samples = out.n_samples;
    best.k_post_mean = out.k_post_mean; best.k_post_sd = out.k_post_sd;
    best.k_eff = 1;

    for (int KK = 1; KK <= kmax; ++KK) {
        double w[kKMax], mu[kKMax], nu[kKMax], k1[kKMax], k2[kKMax], lam[kKMax], lc[kKMax];
        for (int k = 0; k < KK; ++k) {
            if (k < int(seeds.size())) {
                // 用 MCMC 末态簇做初值(位置/精度都对,避免 EM 收敛到"一个巨宽成分")
                w[k] = std::max(seeds[size_t(k)][0], 1e-4);
                mu[k] = seeds[size_t(k)][1]; nu[k] = seeds[size_t(k)][2];
                k1[k] = std::min(std::max(seeds[size_t(k)][3], 0.05), 800.0);
                k2[k] = std::min(std::max(seeds[size_t(k)][4], 0.05), 800.0);
                lam[k] = seeds[size_t(k)][5];
                double cap = 0.95 * std::sqrt(k1[k] * k2[k]);
                if (std::fabs(lam[k]) > cap) lam[k] = (lam[k] >= 0 ? cap : -cap);
                lc[k] = logC_tab(tab, k1[k], k2[k], lam[k]);
                continue;
            }
            double q = (double(k) + 0.5) / double(KK);
            double acc = 0.0; int pick = 0;
            for (int c = 0; c < kCells; ++c) { acc += mass[c]; if (acc >= q) { pick = c; break; } }
            mu[k] = bin_center(pick / kBins);
            nu[k] = bin_center(pick % kBins);
            k1[k] = 4; k2[k] = 4; lam[k] = 0;
            w[k] = 1.0 / double(KK);
            lc[k] = logC_tab(tab, k1[k], k2[k], lam[k]);
        }
        for (int it = 0; it < em_iters; ++it) {
            double den[kCells];
            for (int c = 0; c < kCells; ++c) {
                int ix = c / kBins, iy = c % kBins;
                double lds[16]; double mx = -3e38;
                for (int k = 0; k < KK; ++k) {
                    double cm = std::cos(mu[k]), sm = std::sin(mu[k]);
                    double cn = std::cos(nu[k]), sn = std::sin(nu[k]);
                    double c1 = cph[ix] * cm + sph[ix] * sm, q1 = sph[ix] * cm - cph[ix] * sm;
                    double c2 = cph[iy] * cn + sph[iy] * sn, q2 = sph[iy] * cn - cph[iy] * sn;
                    lds[k] = std::log(std::max(w[k], 1e-30)) + lc[k] + k1[k] * c1 + k2[k] * c2 + lam[k] * q1 * q2;
                    if (lds[k] > mx) mx = lds[k];
                }
                double sum = 0.0;
                for (int k = 0; k < KK; ++k) sum += std::exp(lds[k] - mx);
                den[c] = std::max(sum * std::exp(mx), 1e-300);
            }
            for (int k = 0; k < KK; ++k) {
                double sr = 0, sA = 0, sB = 0, sC = 0, sD = 0, sS = 0, sV1 = 0, sV2 = 0;
                for (int c = 0; c < kCells; ++c) {
                    int ix = c / kBins, iy = c % kBins;
                    double cm = std::cos(mu[k]), sm = std::sin(mu[k]);
                    double cn = std::cos(nu[k]), sn = std::sin(nu[k]);
                    double c1 = cph[ix] * cm + sph[ix] * sm, q1 = sph[ix] * cm - cph[ix] * sm;
                    double c2 = cph[iy] * cn + sph[iy] * sn, q2 = sph[iy] * cn - cph[iy] * sn;
                    double ld = std::log(std::max(w[k], 1e-30)) + lc[k] + k1[k] * c1 + k2[k] * c2 + lam[k] * q1 * q2;
                    double r = std::exp(ld) / den[c] * mass[c];
                    sr += r;
                    sA += r * cph[ix]; sB += r * sph[ix];
                    sC += r * cph[iy]; sD += r * sph[iy];
                    sS += r * q1 * q2;
                    sV1 += r * std::atan2(q1, c1) * std::atan2(q1, c1);
                    sV2 += r * std::atan2(q2, c2) * std::atan2(q2, c2);
                }
                double tr = std::max(sr, 1e-12);
                w[k] = tr;
                double R1 = std::sqrt(sA * sA + sB * sB) / tr;
                double R2 = std::sqrt(sC * sC + sD * sD) / tr;
                // κ 用二阶矩估计(带格内方差修正 h^2/12):合向量法在峰比格子窄时被严重低估,
                // 而 10° 网格上恰好全是窄峰。宽峰时两者一致,取较大者更稳。
                const double hcell = kTwoPi / double(kBins);
                double var1 = sV1 / tr + hcell * hcell / 12.0;
                double var2 = sV2 / tr + hcell * hcell / 12.0;
                double k1v = 1.0 / std::max(var1, 1e-6);
                double k2v = 1.0 / std::max(var2, 1e-6);
                double k1m = std::max(invA(R1), k1v);
                double k2m = std::max(invA(R2), k2v);
                k1m = std::min(k1m, 200.0); k2m = std::min(k2m, 200.0);
                double rho = sS / tr / std::max(std::sqrt(std::max(1 - R1 * R1, 1e-6)) *
                                                std::sqrt(std::max(1 - R2 * R2, 1e-6)), 1e-6);
                rho = std::min(std::max(rho, -0.95), 0.95);
                double lm = rho * std::sqrt(std::max(k1m * k2m, 1e-6));
                mu[k] = std::atan2(sB, sA);
                nu[k] = std::atan2(sD, sC);
                k1[k] = k1m + lm * lm / std::max(k2m, 1e-3);
                k2[k] = k2m + lm * lm / std::max(k1m, 1e-3);
                // 36x36 网格分辨率下 κ 没有意义地大;夹到可行域,避免归一化常数失真
                k1[k] = std::min(std::max(k1[k], 1e-2), 200.0);
                k2[k] = std::min(std::max(k2[k], 1e-2), 200.0);
                double cap = 0.98 * std::sqrt(k1[k] * k2[k]);
                if (std::fabs(lm) > cap) lm = (lm >= 0 ? cap : -cap);
                if (lm * lm > 0.98 * k1[k] * k2[k]) lm = 0.9899 * std::sqrt(k1[k] * k2[k]) * (rho >= 0 ? 1 : -1);
                lam[k] = lm;
                lc[k] = logC_tab(tab, k1[k], k2[k], lam[k]);
            }
        }
        double ll = 0.0;
        for (int c = 0; c < kCells; ++c) {
            int ix = c / kBins, iy = c % kBins;
            double lds[16]; double mx = -3e38;
            for (int k = 0; k < KK; ++k) {
                double cm = std::cos(mu[k]), sm = std::sin(mu[k]);
                double cn = std::cos(nu[k]), sn = std::sin(nu[k]);
                double c1 = cph[ix] * cm + sph[ix] * sm, q1 = sph[ix] * cm - cph[ix] * sm;
                double c2 = cph[iy] * cn + sph[iy] * sn, q2 = sph[iy] * cn - cph[iy] * sn;
                lds[k] = std::log(std::max(w[k], 1e-30)) + lc[k] + k1[k] * c1 + k2[k] * c2 + lam[k] * q1 * q2;
                if (lds[k] > mx) mx = lds[k];
            }
            double sum = 0.0;
            for (int k = 0; k < KK; ++k) sum += std::exp(lds[k] - mx);
            ll += mass[c] * (mx + std::log(std::max(sum, 1e-300)));
        }
        // 网格上的平均对数密度要乘回样本数,才是真正的对数似然
        // (漏了这个 n,BIC 会永远选 K'=1 —— 这就是之前所有切片都被压成单成分的原因)
        double bic = -2.0 * ll * double(std::max(n_obs, 1)) + double(KK) * 6.0 * std::log(std::max(double(n_obs), 2.0));
        if (bic < best_bic) {
            best_bic = bic;
            best.k_eff = uint32_t(KK);
            best.loglik_mean = float(ll);
            best.reserved_f = float(bic);
            // 权重归一
            double ws = 0.0;
            for (int k = 0; k < KK; ++k) ws += w[k];
            if (ws <= 0) ws = 1.0;
            for (int k = 0; k < kKMax; ++k) {
                float* d = best.comp + k * 8;
                if (k < KK) {
                    d[0] = float(w[k] / ws); d[1] = float(mu[k]); d[2] = float(nu[k]);
                    d[3] = float(k1[k]); d[4] = float(k2[k]); d[5] = float(lam[k]);
                    d[6] = float(lc[k]);
                } else {
                    d[0] = 0.f; d[1] = 0.f; d[2] = 0.f; d[3] = 1.f; d[4] = 1.f; d[5] = 0.f; d[6] = 0.f;
                }
                d[7] = 0.f;
            }
        }
    }
    out = best;
}

}  // namespace

// ---------------------------------------------------------------- 对外
bool cpu_fit(const Config& cfg, const PointSet& ps, const CtxStat& st,
             std::vector<ParamRec>& recs, std::vector<float>& pred_grid, std::string& err) {
    static BesselTable tab;
    static bool tab_ready = false;
    if (!tab_ready) { build_bessel_table(tab, 1024, 1e-3, 6000.0, 220); tab_ready = true; }

    FitIn F;
    F.ps = &ps; F.st = &st; F.cfg = &cfg; F.tab = &tab;

    std::vector<CtxItem> items;
    for (int c = 0; c < kNctx; ++c) {
        int n = int(st.n[size_t(c)]);
        if (n == 0) continue;
        int tier = (n >= cfg.tier_a_n) ? 3 : (n >= cfg.tier_b_n ? 2 : 1);
        items.push_back({c, tier, n});
    }
    std::sort(items.begin(), items.end(), [](const CtxItem& a, const CtxItem& b) {
        if (a.tier != b.tier) return a.tier > b.tier;
        return a.n > b.n;
    });
    if (cfg.max_ctx > 0 && int(items.size()) > cfg.max_ctx) items.resize(size_t(cfg.max_ctx));

    recs.assign(size_t(kNctx), ParamRec{});
    pred_grid.assign(size_t(kNctx) * kCells, 0.f);
    for (int c = 0; c < kNctx; ++c) {
        recs[size_t(c)].tier = 0;
        recs[size_t(c)].n_obs = st.n[size_t(c)];
        recs[size_t(c)].k_eff = 0;
        for (int k = 0; k < kKMax; ++k) {
            float* d = recs[size_t(c)].comp + k * 8;
            d[0] = 0.f; d[1] = 0.f; d[2] = 0.f; d[3] = 1.f; d[4] = 1.f; d[5] = 0.f; d[6] = 0.f; d[7] = 0.f;
        }
        if (st.n[size_t(c)] == 0) {
            float* d = recs[size_t(c)].comp;
            d[0] = 1.f; d[1] = 0.f; d[2] = 0.f;
            d[3] = cfg.k10; d[4] = cfg.k20; d[5] = cfg.lam0;
            d[6] = float(logC_exact(cfg.k10, cfg.k20, cfg.lam0));
            recs[size_t(c)].k_eff = 1;
            double lc = logC_exact(cfg.k10, cfg.k20, cfg.lam0);
            for (int c2 = 0; c2 < kCells; ++c2) {
                double ph = bin_center(c2 / kBins), qs = bin_center(c2 % kBins);
                pred_grid[size_t(c) * kCells + c2] =
                    float(std::exp(lc + cfg.k10 * std::cos(ph) + cfg.k20 * std::cos(qs)
                                   + cfg.lam0 * std::sin(ph) * std::sin(qs)));
            }
        }
    }

    std::atomic<size_t> next(0);
    std::atomic<int> done(0);
    int nthreads = std::max(1u, std::thread::hardware_concurrency());
    if (nthreads > 6) nthreads = 6;
    const size_t total = items.size();

    auto worker = [&]() {
        for (;;) {
            size_t i = next.fetch_add(1);
            if (i >= total) break;
            const CtxItem& it = items[i];
            std::vector<float> grid;
            ParamRec rec{};
            std::vector<std::array<double, 6>> seeds;
            run_one_ctx(F, it.ctx, it.tier, grid, rec, seeds);
            rec.k_eff = 0;
            compress_ctx(tab, it.n, grid, cfg.k_max_comp, 40, rec, seeds);
            recs[size_t(it.ctx)] = rec;
            for (int c = 0; c < kCells; ++c) pred_grid[size_t(it.ctx) * kCells + c] = grid[size_t(c)];
            int d = done.fetch_add(1) + 1;
            if (d % 200 == 0) {
                std::printf("\r  MCMC %d/%zu 上下文", d, total);
                std::fflush(stdout);
            }
        }
    };
    std::vector<std::thread> th;
    for (int i = 0; i < nthreads; ++i) th.emplace_back(worker);
    for (auto& t : th) t.join();
    std::printf("\n");

    // 归一化预测网格
    for (int c = 0; c < kNctx; ++c) {
        float* s = pred_grid.data() + size_t(c) * kCells;
        double tot = 0.0;
        for (int i = 0; i < kCells; ++i) { s[i] = float(double(s[i]) * kCellArea); tot += s[i]; }
        float inv = tot > 0 ? float(1.0 / tot) : 0.f;
        for (int i = 0; i < kCells; ++i) s[i] *= inv;
    }
    return true;
}

bool cpu_sample(const Config& cfg, const std::vector<ParamRec>& recs,
                std::vector<uint32_t>& counts, std::string& err) {
    counts.assign(size_t(kNctx) * kCells, 0u);
    std::atomic<int> next(0);
    int nthreads = std::max(1u, std::thread::hardware_concurrency());
    if (nthreads > 6) nthreads = 6;
    auto worker = [&]() {
        for (;;) {
            int ctx = next.fetch_add(1);
            if (ctx >= kNctx) break;
            const ParamRec& R = recs[size_t(ctx)];
            int K = int(R.k_eff);
            if (K < 1) K = 1;
            if (K > kKMax) K = kKMax;
            double cw[kKMax];
            double ws = 0.0;
            for (int k = 0; k < K; ++k) { cw[k] = std::max(double(R.comp[k * 8 + 0]), 0.0); ws += cw[k]; }
            if (ws <= 0) { cw[0] = 1.0; ws = 1.0; }
            for (int k = 0; k < K; ++k) cw[k] /= ws;

            std::vector<uint32_t> local(size_t(kCells), 0u);
            const uint64_t M = cfg.m_samples;
            for (uint64_t s = 0; s < M; ++s) {
                Rng r = make_rng(cfg.seed ^ 0xABCDEFull, uint32_t(ctx), 0u, uint32_t(s & 0xFFFFFFFFu), 0x5EEDu);
                double u = u01(r);
                int pick = K - 1;
                double acc = 0.0;
                for (int k = 0; k < K; ++k) { acc += cw[k]; if (u <= acc) { pick = k; break; } }
                double mu = R.comp[pick * 8 + 1], nu = R.comp[pick * 8 + 2];
                double k1 = R.comp[pick * 8 + 3], k2 = R.comp[pick * 8 + 4], lam = R.comp[pick * 8 + 5];
                double ph = mu, ps = nu;
                for (int it = 0; it < 6; ++it) {
                    double s2 = std::sin(ps - nu);
                    double b1 = lam * s2;
                    ph = sample_vm(mu + std::atan2(b1, k1), std::sqrt(k1 * k1 + b1 * b1), r);
                    double s1 = std::sin(ph - mu);
                    double b2 = lam * s1;
                    ps = sample_vm(nu + std::atan2(b2, k2), std::sqrt(k2 * k2 + b2 * b2), r);
                }
                double a1 = std::fmod(ph + kPi, kTwoPi); if (a1 < 0) a1 += kTwoPi;
                double a2 = std::fmod(ps + kPi, kTwoPi); if (a2 < 0) a2 += kTwoPi;
                int ix = int(a1 / (kTwoPi / double(kBins)));
                int iy = int(a2 / (kTwoPi / double(kBins)));
                if (ix < 0) ix = 0; if (ix >= kBins) ix = kBins - 1;
                if (iy < 0) iy = 0; if (iy >= kBins) iy = kBins - 1;
                local[size_t(ix * kBins + iy)]++;
            }
            for (int c = 0; c < kCells; ++c) counts[size_t(ctx) * kCells + c] = local[size_t(c)];
        }
    };
    std::vector<std::thread> th;
    for (int i = 0; i < nthreads; ++i) th.emplace_back(worker);
    for (auto& t : th) t.join();
    return true;
}

bool cpu_histogram(const PointSet& ps, std::vector<uint32_t>& counts) {
    counts.assign(size_t(kNctx) * kCells, 0u);
    for (size_t i = 0; i < ps.size(); ++i) {
        counts[size_t(ps.ctx[i]) * kCells + size_t(ps.bx[i]) * kBins + ps.by[i]]++;
    }
    return true;
}

}  // namespace vd









namespace vd {

// 把压缩参数在 36x36 格心上求值并归一 -> 格概率(用于和式(7)网格对拍)
void model_grid(const ParamRec& R, float* out) {
    double acc[kCells];
    double tot = 0.0;
    int K = int(R.k_eff);
    if (K < 1) K = 1;
    if (K > kKMax) K = kKMax;
    for (int c = 0; c < kCells; ++c) {
        double ph = bin_center(c / kBins), ps = bin_center(c % kBins);
        double cph = std::cos(ph), sph = std::sin(ph), cps = std::cos(ps), sps = std::sin(ps);
        double v = 0.0;
        for (int k = 0; k < K; ++k) {
            const float* d = R.comp + k * 8;
            if (d[0] <= 0.f) continue;
            double cm = std::cos(d[1]), sm = std::sin(d[1]);
            double cn = std::cos(d[2]), sn = std::sin(d[2]);
            double c1 = cph * cm + sph * sm, q1 = sph * cm - cph * sm;
            double c2 = cps * cn + sps * sn, q2 = sps * cn - cps * sn;
            v += double(d[0]) * std::exp(double(d[6]) + d[3] * c1 + d[4] * c2 + d[5] * q1 * q2);
        }
        v *= kCellArea;
        acc[c] = v;
        tot += v;
    }
    if (tot <= 0) tot = 1.0;
    for (int c = 0; c < kCells; ++c) out[c] = float(acc[c] / tot);
}


bool cpu_recompress(const Config& cfg, std::vector<ParamRec>& recs, const float* pred,
                    double* mean_tv, std::string& err) {
    (void)err;
    static BesselTable tab;
    static bool ready = false;
    if (!ready) { build_bessel_table(tab, 1024, 1e-3, 6000.0, 220); ready = true; }

    double sum_tv = 0.0;
    int cnt = 0;
    std::vector<float> mg;
    mg.resize(size_t(kCells));
    for (int c = 0; c < kNctx; ++c) {
        ParamRec src = recs[size_t(c)];
        if (src.n_obs == 0) continue;
        std::vector<float> grid(pred + size_t(c) * kCells, pred + size_t(c + 1) * kCells);
        ParamRec out = src;
        out.k_eff = 0;
        std::vector<std::array<double, 6>> noSeed;
        compress_ctx(tab, int(src.n_obs), grid, cfg.k_max_comp, 60, out, noSeed);
        model_grid(out, mg.data());
        double tv = 0.0;
        double s1 = 0.0, s2 = 0.0;
        for (int i = 0; i < kCells; ++i) { s1 += grid[size_t(i)]; s2 += mg[size_t(i)]; }
        for (int i = 0; i < kCells; ++i) {
            double a = grid[size_t(i)] / std::max(s1, 1e-30);
            double b = mg[size_t(i)] / std::max(s2, 1e-30);
            tv += std::fabs(a - b);
        }
        sum_tv += 0.5 * tv;
        ++cnt;
        recs[size_t(c)] = out;
        if (cnt % 400 == 0) { std::printf("\r  压缩 %d", cnt); std::fflush(stdout); }
    }
    std::printf("\n");
    if (mean_tv) *mean_tv = cnt ? sum_tv / cnt : 0.0;
    return true;
}

}  // namespace vd



