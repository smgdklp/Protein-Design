// von_dpmm / gpu.cu — CUDA 核:直方图、DPMM-MCMC、后验预测网格、压缩、离散采样
#include "gpu.h"

#include <cstdio>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <vector>
#include <cuda_runtime.h>

namespace vd {

// ================================================================ 设备状态
namespace {

struct Dev {
    float*  tab = nullptr;      // Bessel 表
    int     nk = 0, mmax = 0;
    float   logkmin = 0.f, dlogk = 1.f;

    uint32_t* d_start = nullptr;
    uint32_t* d_order = nullptr;
    float *d_cphi = nullptr, *d_sphi = nullptr, *d_cpsi = nullptr, *d_spsi = nullptr;
    uint8_t *d_bx = nullptr, *d_by = nullptr;
    int*     d_label = nullptr;

    uint32_t* d_hist = nullptr;     // [kNctx*kCells]
    ParamRec* d_recs = nullptr;     // [kNctx]
    float*    d_grid = nullptr;     // [kNctx*kCells]

    size_t cap_points = 0;
    size_t free_before = 0, total = 0;
    bool   ready = false;
};
Dev g;

__constant__ float c_lbinom[256];   // log C(2m,m), m<256

void build_lbinom() {
    std::vector<float> v(256);
    for (int m = 0; m < 256; ++m) {
        double s = std::lgamma(2.0 * m + 1.0) - 2.0 * std::lgamma(double(m) + 1.0);
        v[m] = float(s);
    }
    cudaMemcpyToSymbol(c_lbinom, v.data(), v.size() * sizeof(float));
}

inline const char* estr(cudaError_t e) {
    const char* s = cudaGetErrorString(e);
    return s ? s : "unknown-cuda-error";
}
inline void chk(const char* what) {
    cudaError_t e = cudaGetLastError();
    if (e != cudaSuccess) {
        std::fprintf(stderr, "[cuda] %s: %s\n", what, estr(e));
    }
}

}  // namespace

// ================================================================ 设备数学
namespace {

struct Rng {
    uint64_t s;
};

__device__ __forceinline__ uint64_t mix64(uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}
__device__ __forceinline__ Rng make_rng(uint64_t seed, uint32_t ctx, uint32_t chain,
                                        uint32_t iter, uint32_t salt) {
    Rng r;
    uint64_t s = seed
               ^ (uint64_t(ctx + 1) * 0x9E3779B97F4A7C15ULL)
               ^ (uint64_t(chain + 1) * 0xD1B54A32D192ED03ULL)
               ^ (uint64_t(iter + 1) * 0x94D049BB133111EBULL)
               ^ (uint64_t(salt + 1) * 0xA24BAED4963EE407ULL);
    r.s = mix64(s);
    return r;
}
__device__ __forceinline__ float u01(Rng& r) {
    r.s = mix64(r.s);
    return float((r.s >> 11) & 0x1FFFFFFFFFFFFFULL) * (1.0f / 9007199254740992.0f);
}
__device__ __forceinline__ float nrm(Rng& r) {
    float u1 = fmaxf(u01(r), 1e-9f), u2 = u01(r);
    return sqrtf(-2.0f * __logf(u1)) * __cosf(6.28318530718f * u2);
}

// log(e^-k I_m(k)),线性插值查表
__device__ __forceinline__ float logIe(const float* tab, int nk, int mmax,
                                       float logkmin, float dlogk, int m, float k) {
    float lk = __logf(fmaxf(k, 1e-9f));
    float t = (lk - logkmin) / dlogk;
    if (t < 0.f) t = 0.f;
    float tm = float(nk - 1);
    if (t > tm) t = tm;
    int i0 = int(t);
    int i1 = i0 + 1 < nk ? i0 + 1 : nk - 1;
    float f = t - float(i0);
    const float* r0 = tab + size_t(i0) * size_t(mmax + 1) + m;
    const float* r1 = tab + size_t(i1) * size_t(mmax + 1) + m;
    return r0[0] * (1.f - f) + r1[0] * f;
}

// vM sine 的 log 归一化常数
__device__ __forceinline__ float logC(const float* tab, int nk, int mmax,
                                      float logkmin, float dlogk,
                                      float k1, float k2, float lam) {
    k1 = fmaxf(k1, 1e-7f);
    k2 = fmaxf(k2, 1e-7f);
    float t = lam * lam / (4.f * k1 * k2);
    if (t > 0.2499999f) t = 0.2499999f;
    float lt = (t > 1e-30f) ? __logf(t) : -1e30f;
    int M = int(8.f * sqrtf(fmaxf(k1, k2))) + 40;
    if (M > 240) M = 240;
    float best = -3e38f;
    float sum = 0.f;
    for (int m = 0; m <= M; ++m) {
        float s = c_lbinom[m];
        if (m > 0) {
            if (lt < -1e29f) break;
            s += float(m) * lt;
        }
        s += logIe(tab, nk, mmax, logkmin, dlogk, m, k1);
        s += logIe(tab, nk, mmax, logkmin, dlogk, m, k2);
        if (s > best) { sum = sum * __expf(best - s) + 1.f; best = s; }
        else sum += __expf(s - best);
    }
    return -1.3862943611198906f /*log(4pi^2)=log(39.478)*/ - (k1 + k2) - (best + __logf(sum));
}

// 单变量 von Mises 采样(Best & Fisher 1979)
__device__ __forceinline__ float sample_vm(float mu, float kappa, Rng& r) {
    kappa = fabsf(kappa);
    if (!(kappa > 1e-6f)) return mu + 6.28318530718f * u01(r);
    float r_ = 1.f + sqrtf(1.f + 4.f * kappa * kappa);
    float rho = (r_ - sqrtf(2.f * r_)) / (2.f * kappa);
    for (int guard = 0; guard < 64; ++guard) {
        float u1 = u01(r), u2 = u01(r), u3 = u01(r);
        float z = __cosf(3.14159265358979f * u1);
        float f = (1.f + rho * z) / (rho + z);
        float c = kappa * (rho - f);
        bool ok = (c * (2.f - c) - u2 > 0.f);
        if (!ok && c > 0.f) ok = (__logf(c / fmaxf(u2, 1e-12f)) + 1.f - c >= 0.f);
        if (ok) {
            float sg = (u3 - 0.5f >= 0.f) ? 1.f : -1.f;
            float a = fminf(fmaxf(f, -1.f), 1.f);
            return mu + sg * acosf(a);
        }
    }
    return mu;
}

// 逆 Bessel 比:由合向量长度 R 估计 kappa
__device__ __forceinline__ float invA(float R) {
    R = fminf(fmaxf(R, 0.f), 0.999999f);
    if (R < 0.53f) return 2.f * R + R * R * R + 5.f * R * R * R * R * R / 6.f;
    if (R < 0.85f) return -0.4f + 1.39f * R + 0.43f / (1.f - R);
    return 1.f / fmaxf(R * R * R - 4.f * R * R + 3.f * R, 1e-6f);
}

__device__ __forceinline__ void bartlett_wishart(Rng& r, float nu, float b11, float b21, float b22,
                                                 float& o11, float& o12, float& o22) {
    // 2x2 Wishart(nu, B): Omega = L L^T,L 为下三角
    // L11 = sqrt(B11 * chi2_nu), L21 = B21*chi_nu/sqrt(chi2_nu) + sqrt(B22 - B21^2/B11)*N(0,1)
    // 简化:chi2_nu = sum of nu squared normals (nu 取整)
    int n = int(nu + 0.5f);
    float c2 = 0.f;
    for (int i = 0; i < n; ++i) { float z = nrm(r); c2 += z * z; }
    if (c2 < 1e-9f) c2 = 1e-9f;
    float l11 = sqrtf(fmaxf(b11 * c2, 1e-12f));
    float z2 = nrm(r);
    float cond = fmaxf(b22 - b21 * b21 / fmaxf(b11, 1e-9f), 1e-9f);
    float l21 = b21 * sqrtf(c2) / l11 + sqrtf(cond) * z2;
    float l22 = sqrtf(fmaxf(cond * (c2 - z2 * z2) , 1e-6f));
    (void)l22;
    o11 = l11 * l11;
    o12 = l11 * l21;
    o22 = l21 * l21 + l22 * l22;
}

// 块内规约(bx: 至少 32 float 的共享缓冲)
__device__ __forceinline__ float block_sum(float v, float* shm) {
    int lane = threadIdx.x & 31, warp = threadIdx.x >> 5;
    for (int o = 16; o > 0; o >>= 1) v += __shfl_down_sync(0xffffffffu, v, o);
    if (lane == 0) shm[warp] = v;
    __syncthreads();
    int nw = (blockDim.x + 31) >> 5;
    v = (threadIdx.x < nw) ? shm[threadIdx.x] : 0.f;
    if (warp == 0) {
        for (int o = 16; o > 0; o >>= 1) v += __shfl_down_sync(0xffffffffu, v, o);
    }
    __syncthreads();
    return v;
}

__device__ __forceinline__ float warp_sum(float v) {
    for (int o = 16; o > 0; o >>= 1) v += __shfl_down_sync(0xffffffffu, v, o);
    return v;
}

}  // namespace

// ================================================================ 直方图
namespace {
__global__ void k_hist(const uint8_t* __restrict__ bx, const uint8_t* __restrict__ by,
                       const uint32_t* __restrict__ ctx, size_t n, uint32_t* __restrict__ out) {
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    size_t stride = size_t(gridDim.x) * blockDim.x;
    for (; i < n; i += stride) {
        uint32_t cell = bx[i] * uint32_t(kBins) + by[i];
        atomicAdd(&out[size_t(ctx[i]) * kCells + cell], 1u);
    }
}
}  // namespace

// ================================================================ MCMC
namespace {

struct FitArgs {
    const uint32_t* start;
    const uint32_t* order;
    const float *cphi, *sphi, *cpsi, *spsi;
    int* label;
    const float* tab;
    int nk, mmax;
    float logkmin, dlogk;
    const uint32_t* ctx_list;
    const uint8_t* tier_list;
    float tau0, k10, k20, lam0;
    int nu0;
    float b0;
    int k_aux;
    uint64_t seed;
    int iters_a, burn_a, thin_a, chains_a;
    int iters_b, burn_b, thin_b, chains_b;
    int iters_c, burn_c, thin_c, chains_c;
    ParamRec* recs;
    float* grid;
};

__device__ __forceinline__ void comp_logdens(const float* tab, int nk, int mmax,
                                             float logkmin, float dlogk,
                                             float cph, float sph, float cps, float sps,
                                             float mu, float nu, float k1, float k2,
                                             float lam, float lc, float& out) {
    float cm = __cosf(mu), sm = __sinf(mu);
    float cn = __cosf(nu), sn = __sinf(nu);
    float c1 = cph * cm + sph * sm;          // cos(phi-mu)
    float s1 = sph * cm - cph * sm;          // sin(phi-mu)
    float c2 = cps * cn + sps * sn;          // cos(psi-nu)
    float s2 = sps * cn - cps * sn;          // sin(psi-nu)
    out = lc + k1 * c1 + k2 * c2 + lam * s1 * s2;
    (void)tab; (void)nk; (void)mmax; (void)logkmin; (void)dlogk;
}

__global__ void k_mcmc(FitArgs A) {
    const int BS = blockDim.x;
    const int tid = threadIdx.x;
    const int lane = tid & 31, warp = tid >> 5;
    const int nwarps = (BS + 31) >> 5;

    const uint32_t ctx = A.ctx_list[blockIdx.x];
    const int tier = int(A.tier_list[blockIdx.x]);
    const uint32_t p0 = A.start[ctx];
    const int n = int(A.start[ctx + 1] - p0);
    int* labp = A.label + p0;

    int iters, burn, thin, chains;
    if (tier >= 3)      { iters = A.iters_a; burn = A.burn_a; thin = A.thin_a; chains = A.chains_a; }
    else if (tier == 2) { iters = A.iters_b; burn = A.burn_b; thin = A.thin_b; chains = A.chains_b; }
    else                { iters = A.iters_c; burn = A.burn_c; thin = A.thin_c; chains = A.chains_c; }

    __shared__ float s_mu[kKCap], s_nu[kKCap], s_k1[kKCap], s_k2[kKCap], s_lam[kKCap], s_lc[kKCap];
    __shared__ int   s_cnt[kKCap];
    __shared__ float s_acc[kKCap * 3];
    __shared__ float s_amu[8], s_anu[8], s_ak1[8], s_ak2[8], s_alam[8], s_alc[8];
    __shared__ float s_grid[kCells];
    __shared__ float s_red[32];
    __shared__ float s_scale[3];
    __shared__ int   s_K, s_Kaux;
    __shared__ float s_dummy;

    const float K10 = A.k10, K20 = A.k20, LAM0 = A.lam0;
    const int KAUX = A.k_aux > 8 ? 8 : A.k_aux;

    // 初始化共享
    for (int c = tid; c < kCells; c += BS) s_grid[c] = 0.f;
    if (tid == 0) { s_scale[0] = 0.35f; s_scale[1] = 0.35f; s_scale[2] = 0.35f; }
    __syncthreads();

    float km_sum = 0.f, km2_sum = 0.f;
    int   nkeep = 0;

    for (int chain = 0; chain < chains; ++chain) {
        // ---------- 初始化:全一簇
        if (tid == 0) {
            s_K = 1;
            s_cnt[0] = n;
            s_mu[0] = 0.f; s_nu[0] = 0.f;
            s_k1[0] = 1.f; s_k2[0] = 1.f; s_lam[0] = 0.f;
            s_lc[0] = logC(A.tab, A.nk, A.mmax, A.logkmin, A.dlogk, 1.f, 1.f, 0.f);
            // 用矩估计给一个像样的起点
        }
        __syncthreads();
        for (int t = tid; t < n; t += BS) labp[t] = 0;
        __syncthreads();

        for (int it = 0; it < iters; ++it) {
            const int K = s_K;
            // ---------- 抽辅助簇
            if (tid < KAUX) {
                Rng r = make_rng(A.seed, ctx, uint32_t(chain), uint32_t(it), uint32_t(tid) + 17u);
                float mu = sample_vm(0.f, K10, r);
                float nu = sample_vm(0.f, K20, r);
                float o11, o12, o22;
                bartlett_wishart(r, float(A.nu0), A.b0, 0.f, A.b0, o11, o12, o22);
                float k1 = o11, k2 = fmaxf(o22, 1e-4f), lam = -o12;
                if (lam * lam > 0.98f * k1 * k2) lam = sqrtf(0.98f * k1 * k2) * (lam >= 0 ? 1.f : -1.f);
                s_amu[tid] = mu; s_anu[tid] = nu; s_ak1[tid] = k1; s_ak2[tid] = k2; s_alam[tid] = lam;
                s_alc[tid] = logC(A.tab, A.nk, A.mmax, A.logkmin, A.dlogk, k1, k2, lam);
            }
            __syncthreads();

            const float lden = __logf(float(n) + A.tau0);
            const float lw_aux = __logf(A.tau0 / float(KAUX)) - lden;

            // ---------- 分配(Collapsed 之外:给定成分参数直接算似然)
            for (int t = tid; t < n; t += BS) {
                uint32_t gp = p0 + uint32_t(t);
                float cph = A.cphi[gp], sph = A.sphi[gp], cps = A.cpsi[gp], sps = A.spsi[gp];
                float best = -3e38f;
                float sc[32 + 8];
                int ncomp = K + KAUX;
                for (int k = 0; k < ncomp; ++k) {
                    float mu, nu, k1, k2, lam, lc, lw;
                    if (k < K) { mu = s_mu[k]; nu = s_nu[k]; k1 = s_k1[k]; k2 = s_k2[k]; lam = s_lam[k]; lc = s_lc[k];
                                 lw = __logf(float(s_cnt[k])) - lden; }
                    else { int j = k - K; mu = s_amu[j]; nu = s_anu[j]; k1 = s_ak1[j]; k2 = s_ak2[j]; lam = s_alam[j]; lc = s_alc[j]; lw = lw_aux; }
                    float ld;
                    comp_logdens(A.tab, A.nk, A.mmax, A.logkmin, A.dlogk, cph, sph, cps, sps,
                                 mu, nu, k1, k2, lam, lc, ld);
                    ld += lw;
                    sc[k] = ld;
                    if (ld > best) best = ld;
                }
                float sum = 0.f;
                for (int k = 0; k < ncomp; ++k) { sc[k] = __expf(sc[k] - best); sum += sc[k]; }
                Rng r = make_rng(A.seed, ctx, uint32_t(chain), uint32_t(it), uint32_t(t) * 3u + 1013u);
                float u = u01(r) * sum;
                int pick = ncomp - 1;
                float acc = 0.f;
                for (int k = 0; k < ncomp; ++k) { acc += sc[k]; if (u <= acc) { pick = k; break; } }
                labp[t] = pick;
            }
            __syncthreads();

            // ---------- 落簇:辅助簇变成真实簇(块更新)
            if (tid == 0) {
                for (int j = 0; j < KAUX; ++j) {
                    int c = 0;
                    for (int t = 0; t < n; ++t) if (labp[t] == K + j) ++c;
                    if (c > 0) {
                        if (s_K < kKCap) {
                            int k = s_K++;
                            s_mu[k] = s_amu[j]; s_nu[k] = s_anu[j];
                            s_k1[k] = s_ak1[j]; s_k2[k] = s_ak2[j]; s_lam[k] = s_alam[j]; s_lc[k] = s_alc[j];
                            s_cnt[k] = c;
                            for (int t = 0; t < n; ++t) if (labp[t] == K + j) labp[t] = k;
                        } else {
                            for (int t = 0; t < n; ++t) if (labp[t] == K + j) labp[t] = 0;
                        }
                    }
                }
            }
            __syncthreads();

            const int K2 = s_K;
            // ---------- 重数
            for (int k = tid; k < K2; k += BS) {
                int c = 0;
                for (int t = 0; t < n; ++t) if (labp[t] == k) ++c;
                s_cnt[k] = c;
            }
            __syncthreads();

            // ---------- 均值:μ|ν (精确 vM)
            for (int i = 0; i < kKCap * 3; ++i) if (i % kKCap == tid) { s_acc[i] = 0.f; }
            __syncthreads();
            for (int kk = warp; kk < K2; kk += nwarps) {
                float sa = 0.f, sb = 0.f;
                float k1 = s_k1[kk], lam = s_lam[kk], nu = s_nu[kk];
                for (int t = lane; t < n; t += 32) {
                    if (labp[t] != kk) continue;
                    uint32_t gp = p0 + uint32_t(t);
                    float cph = A.cphi[gp], sph = A.sphi[gp];
                    float cps = A.cpsi[gp], sps = A.spsi[gp];
                    float sn = __sinf(nu), cn = __cosf(nu);
                    float s2 = sps * cn - cps * sn;      // sin(psi-nu)
                    sa += k1 * cph + lam * s2 * sph;
                    sb += k1 * sph - lam * s2 * cph;
                }
                sa = warp_sum(sa); sb = warp_sum(sb);
                if (lane == 0) { atomicAdd(&s_acc[kk], sa); atomicAdd(&s_acc[kKCap + kk], sb); }
            }
            __syncthreads();
            if (tid < K2) {
                float Aa = s_acc[tid] + K10;                       // 先验:mu0=0
                float Bb = s_acc[kKCap + tid] - LAM0 * __sinf(s_nu[tid]);
                float kap = sqrtf(Aa * Aa + Bb * Bb);
                float mu0 = atan2f(Bb, Aa);
                Rng r = make_rng(A.seed, ctx, uint32_t(chain), uint32_t(it), uint32_t(tid) + 7001u);
                s_mu[tid] = sample_vm(mu0, kap, r);
            }
            __syncthreads();

            // ---------- 均值:ν|μ
            for (int i = 0; i < kKCap * 3; ++i) if (i % kKCap == tid) { s_acc[i] = 0.f; }
            __syncthreads();
            for (int kk = warp; kk < K2; kk += nwarps) {
                float sa = 0.f, sb = 0.f;
                float k2 = s_k2[kk], lam = s_lam[kk], mu = s_mu[kk];
                for (int t = lane; t < n; t += 32) {
                    if (labp[t] != kk) continue;
                    uint32_t gp = p0 + uint32_t(t);
                    float cph = A.cphi[gp], sph = A.sphi[gp];
                    float cps = A.cpsi[gp], sps = A.spsi[gp];
                    float sm = __sinf(mu), cm = __cosf(mu);
                    float s1 = sph * cm - cph * sm;      // sin(phi-mu)
                    sa += k2 * cps + lam * s1 * sps;
                    sb += k2 * sps - lam * s1 * cps;
                }
                sa = warp_sum(sa); sb = warp_sum(sb);
                if (lane == 0) { atomicAdd(&s_acc[kk], sa); atomicAdd(&s_acc[kKCap + kk], sb); }
            }
            __syncthreads();
            if (tid < K2) {
                float Aa = s_acc[tid] + K20;
                float Bb = s_acc[kKCap + tid] - LAM0 * __sinf(s_mu[tid]);
                float kap = sqrtf(Aa * Aa + Bb * Bb);
                float nu0 = atan2f(Bb, Aa);
                Rng r = make_rng(A.seed, ctx, uint32_t(chain), uint32_t(it), uint32_t(tid) + 9001u);
                s_nu[tid] = sample_vm(nu0, kap, r);
            }
            __syncthreads();

            // ---------- 精度:(κ1,κ2,λ) 随机游走 MH(Wishart 先验)
            for (int i = 0; i < kKCap * 3; ++i) if (i % kKCap == tid) { s_acc[i] = 0.f; }
            __syncthreads();
            for (int kk = warp; kk < K2; kk += nwarps) {
                float s1 = 0.f, s2 = 0.f, s3 = 0.f;
                float mu = s_mu[kk], nu = s_nu[kk];
                float cm = __cosf(mu), sm = __sinf(mu), cn = __cosf(nu), sn = __sinf(nu);
                for (int t = lane; t < n; t += 32) {
                    if (labp[t] != kk) continue;
                    uint32_t gp = p0 + uint32_t(t);
                    float cph = A.cphi[gp], sph = A.sphi[gp];
                    float cps = A.cpsi[gp], sps = A.spsi[gp];
                    float c1 = cph * cm + sph * sm;
                    float q1 = sph * cm - cph * sm;
                    float c2 = cps * cn + sps * sn;
                    float q2 = sps * cn - cps * sn;
                    s1 += c1; s2 += c2; s3 += q1 * q2;
                }
                s1 = warp_sum(s1); s2 = warp_sum(s2); s3 = warp_sum(s3);
                if (lane == 0) {
                    atomicAdd(&s_acc[kk], s1);
                    atomicAdd(&s_acc[kKCap + kk], s2);
                    atomicAdd(&s_acc[2 * kKCap + kk], s3);
                }
            }
            __syncthreads();
            if (tid < K2) {
                float S1 = s_acc[tid], S2 = s_acc[kKCap + tid], S3 = s_acc[2 * kKCap + tid];
                float cnt = float(s_cnt[tid]);
                float k1 = s_k1[tid], k2 = s_k2[tid], lam = s_lam[tid];
                float lc0 = s_lc[tid];
                float lp0 = k1 * S1 + k2 * S2 + lam * S3 + cnt * lc0
                          - 0.5f * __logf(fmaxf(k1 * k2 - lam * lam, 1e-12f))
                          - 2.f * (1.f / A.b0) * 0.5f * (k1 + k2) * A.b0;  // tr(B0^-1 W)/2, B0=b0*I
                Rng r = make_rng(A.seed, ctx, uint32_t(chain), uint32_t(it), uint32_t(tid) + 12007u);
                float rho = lam / sqrtf(fmaxf(k1 * k2, 1e-12f));
                float z = atanhf(fminf(fmaxf(rho, -0.999f), 0.999f));
                float nk1 = k1 * __expf(s_scale[0] * nrm(r));
                float nk2 = k2 * __expf(s_scale[1] * nrm(r));
                float nz = z + s_scale[2] * nrm(r);
                float nrho = tanhf(fminf(fmaxf(nz, -3.f), 3.f));
                float nlam = nrho * sqrtf(nk1 * nk2);
                float nlc = logC(A.tab, A.nk, A.mmax, A.logkmin, A.dlogk, nk1, nk2, nlam);
                float lp1 = nk1 * S1 + nk2 * S2 + nlam * S3 + cnt * nlc
                          - 0.5f * __logf(fmaxf(nk1 * nk2 - nlam * nlam, 1e-12f))
                          - 2.f * (1.f / A.b0) * 0.5f * (nk1 + nk2) * A.b0;
                // Jacobian:k1,e^{x}; 对称提议下 log-k 空间的 Jacobian 抵消
                float u = u01(r);
                if (__logf(fmaxf(u, 1e-12f)) < (lp1 - lp0)) {
                    s_k1[tid] = nk1; s_k2[tid] = nk2; s_lam[tid] = nlam; s_lc[tid] = nlc;
                }
                // 自适应(仅 burn-in 后期)
                if (it < burn && it > burn / 2 && it % 50 == 0 && tid == 0) {
                    // 粗调:接受率无法在此统计,使用固定衰减
                }
                (void)s_dummy;
            }
            __syncthreads();

            // ---------- 网格累积(式 7)
            if (it >= burn && ((it - burn) % thin) == 0) {
                const int Kc = s_K;
                float wsum = 0.f;
                for (int k = 0; k < Kc; ++k) wsum += float(s_cnt[k]);
                float wnew = A.tau0 / (wsum + A.tau0);
                float lnew = logC(A.tab, A.nk, A.mmax, A.logkmin, A.dlogk, K10, K20, LAM0);
                for (int c = tid; c < kCells; c += BS) {
                    int ix = c / kBins, iy = c % kBins;
                    float ph = -3.14159265358979f + (float(ix) + 0.5f) * 6.28318530718f / float(kBins);
                    float ps = -3.14159265358979f + (float(iy) + 0.5f) * 6.28318530718f / float(kBins);
                    float cph = __cosf(ph), sph = __sinf(ph), cps = __cosf(ps), sps = __sinf(ps);
                    float v = 0.f;
                    for (int k = 0; k < Kc; ++k) {
                        float cm = __cosf(s_mu[k]), sm = __sinf(s_mu[k]);
                        float cn = __cosf(s_nu[k]), sn = __sinf(s_nu[k]);
                        float c1 = cph * cm + sph * sm, q1 = sph * cm - cph * sm;
                        float c2 = cps * cn + sps * sn, q2 = sps * cn - cps * sn;
                        float w = float(s_cnt[k]) / (wsum + A.tau0);
                        v += w * __expf(s_lc[k] + s_k1[k] * c1 + s_k2[k] * c2 + s_lam[k] * q1 * q2);
                    }
                    // 新簇(先验预测)项
                    {
                        float cm = 1.f, sm = 0.f, cn = 1.f, sn = 0.f;
                        float c1 = cph * cm + sph * sm, q1 = sph * cm - cph * sm;
                        float c2 = cps * cn + sps * sn, q2 = sps * cn - cps * sn;
                        v += wnew * __expf(lnew + K10 * c1 + K20 * c2 + LAM0 * q1 * q2);
                    }
                    s_grid[c] += v;
                }
                if (tid == 0) { km_sum += float(Kc); km2_sum += float(Kc) * float(Kc); ++nkeep; }
            }
            __syncthreads();
        }
    }

    // ---------- 写回网格(未归一,密度和)
    for (int c = tid; c < kCells; c += BS) A.grid[size_t(ctx) * kCells + c] = s_grid[c];

    ParamRec& R = A.recs[ctx];
    if (tid == 0) {
        R.tier = uint32_t(tier);
        R.n_obs = uint32_t(n);
        R.n_samples = uint32_t(nkeep);
        R.k_eff = uint32_t(s_K);
        float m = nkeep > 0 ? km_sum / float(nkeep) : 0.f;
        float m2 = nkeep > 0 ? km2_sum / float(nkeep) : 0.f;
        R.k_post_mean = m;
        R.k_post_sd = sqrtf(fmaxf(m2 - m * m, 0.f));
        R.loglik_mean = 0.f;
        R.reserved_f = 0.f;
        for (int k = 0; k < kKMax; ++k) {
            float* d = R.comp + k * 8;
            if (k < s_K) {
                d[0] = float(s_cnt[k]) / (float(n) + A.tau0);
                d[1] = s_mu[k]; d[2] = s_nu[k];
                d[3] = s_k1[k]; d[4] = s_k2[k]; d[5] = s_lam[k]; d[6] = s_lc[k];
            } else {
                d[0] = 0.f; d[1] = 0.f; d[2] = 0.f; d[3] = 1.f; d[4] = 1.f; d[5] = 0.f; d[6] = 0.f;
            }
            d[7] = 0.f;
        }
    }
}

// ================================================================ 压缩(网格 -> K' 个 vM sine)
__global__ void k_compress(const float* __restrict__ grid, ParamRec* __restrict__ recs,
                           const float* tab, int nk, int mmax, float logkmin, float dlogk,
                           int kmax, int em_iters) {
    const int ctx = blockIdx.x;
    const int tid = threadIdx.x;
    const int BS = blockDim.x;
    const int nwarps = (BS + 31) >> 5;
    const int lane = tid & 31, warp = tid >> 5;

    __shared__ float s_cph[kBins], s_sph[kBins], s_cps[kBins], s_sps[kBins];
    __shared__ float s_mass[kCells];
    __shared__ float s_den[kCells];
    __shared__ float s_w[kKMax], s_mu[kKMax], s_nu[kKMax], s_k1[kKMax], s_k2[kKMax], s_lam[kKMax], s_lc[kKMax];
    __shared__ float s_red[32];

    const float* G = grid + size_t(ctx) * kCells;
    float tot = 0.f;
    for (int c = tid; c < kCells; c += BS) tot += G[c];
    tot = block_sum(tot, s_red);
    if (tot <= 0.f) tot = 1.f;
    for (int c = tid; c < kCells; c += BS) s_mass[c] = G[c] / tot;

    for (int k = tid; k < kBins; k += BS) {
        float ph = -3.14159265358979f + (float(k) + 0.5f) * 6.28318530718f / float(kBins);
        s_cph[k] = __cosf(ph); s_sph[k] = __sinf(ph);
        s_cps[k] = s_cph[k];    s_sps[k] = s_sph[k];
    }
    __syncthreads();

    float best_bic = 3e38f;
    int   best_k = 1;
    // 保存最优结果
    ParamRec bestrec;

    for (int KK = 1; KK <= kmax; ++KK) {
        // ---- 初始化:KMeans++ 风格(确定性)
        if (tid == 0) {
            for (int k = 0; k < KK; ++k) {
                Rng r = make_rng(777u, uint32_t(ctx), 0u, uint32_t(k), 4242u);
                // 用累积质量分位点初始化
                float q = (float(k) + 0.5f) / float(KK);
                float acc = 0.f;
                int pick = 0;
                for (int c = 0; c < kCells; ++c) { acc += s_mass[c]; if (acc >= q) { pick = c; break; } }
                int ix = pick / kBins, iy = pick % kBins;
                s_mu[k] = -3.14159265358979f + (float(ix) + 0.5f) * 6.28318530718f / float(kBins);
                s_nu[k] = -3.14159265358979f + (float(iy) + 0.5f) * 6.28318530718f / float(kBins);
                s_k1[k] = 4.f; s_k2[k] = 4.f; s_lam[k] = 0.f;
                s_w[k] = 1.f / float(KK);
                s_lc[k] = logC(tab, nk, mmax, logkmin, dlogk, s_k1[k], s_k2[k], s_lam[k]);
                (void)r;
            }
        }
        __syncthreads();

        for (int it = 0; it < em_iters; ++it) {
            // E 步分母
            for (int c = tid; c < kCells; c += BS) {
                int ix = c / kBins, iy = c % kBins;
                float cph = s_cph[ix], sph = s_sph[ix], cps = s_cps[iy], sps = s_sps[iy];
                float mx = -3e38f;
                float tmp[16];
                for (int k = 0; k < KK; ++k) {
                    float cm = __cosf(s_mu[k]), sm = __sinf(s_mu[k]);
                    float cn = __cosf(s_nu[k]), sn = __sinf(s_nu[k]);
                    float c1 = cph * cm + sph * sm, q1 = sph * cm - cph * sm;
                    float c2 = cps * cn + sps * sn, q2 = sps * cn - cps * sn;
                    tmp[k] = __logf(fmaxf(s_w[k], 1e-30f)) + s_lc[k] + s_k1[k] * c1 + s_k2[k] * c2 + s_lam[k] * q1 * q2;
                    if (tmp[k] > mx) mx = tmp[k];
                }
                float sum = 0.f;
                for (int k = 0; k < KK; ++k) { tmp[k] = __expf(tmp[k] - mx); sum += tmp[k]; }
                s_den[c] = sum * __expf(mx);
            }
            __syncthreads();

            // M 步:逐成分统计
            for (int k = 0; k < KK; ++k) {
                float sr = 0.f, sA = 0.f, sB = 0.f, sC = 0.f, sD = 0.f, sS = 0.f;
                for (int c = tid; c < kCells; c += BS) {
                    int ix = c / kBins, iy = c % kBins;
                    float cph = s_cph[ix], sph = s_sph[ix], cps = s_cps[iy], sps = s_sps[iy];
                    float cm = __cosf(s_mu[k]), sm = __sinf(s_mu[k]);
                    float cn = __cosf(s_nu[k]), sn = __sinf(s_nu[k]);
                    float c1 = cph * cm + sph * sm, q1 = sph * cm - cph * sm;
                    float c2 = cps * cn + sps * sn, q2 = sps * cn - cps * sn;
                    float ld = __logf(fmaxf(s_w[k], 1e-30f)) + s_lc[k] + s_k1[k] * c1 + s_k2[k] * c2 + s_lam[k] * q1 * q2;
                    float r = __expf(ld) / fmaxf(s_den[c], 1e-30f) * s_mass[c];
                    sr += r;
                    sA += r * cph; sB += r * sph;
                    sC += r * cps; sD += r * sps;
                    sS += r * q1 * q2;
                }
                float v[6] = {sr, sA, sB, sC, sD, sS};
                float sh[6][32];
                (void)sh;
                for (int j = 0; j < 6; ++j) {
                    float x = block_sum(v[j], s_red);
                    v[j] = x;
                }
                if (tid == 0) {
                    float tot_r = fmaxf(v[0], 1e-12f);
                    s_w[k] = tot_r;
                    float mu = atan2f(v[2], v[1]);
                    float nu = atan2f(v[4], v[3]);
                    float R1 = sqrtf(v[1] * v[1] + v[2] * v[2]) / tot_r;
                    float R2 = sqrtf(v[3] * v[3] + v[4] * v[4]) / tot_r;
                    float k1m = invA(R1), k2m = invA(R2);
                    float rho = v[5] / tot_r / fmaxf(sqrtf(fmaxf(1.f - R1 * R1, 1e-6f)) *
                                                     sqrtf(fmaxf(1.f - R2 * R2, 1e-6f)), 1e-6f);
                    rho = fminf(fmaxf(rho, -0.95f), 0.95f);
                    float lam = rho * sqrtf(fmaxf(k1m * k2m, 1e-6f));
                    float k1 = k1m + lam * lam / fmaxf(k2m, 1e-3f);
                    float k2 = k2m + lam * lam / fmaxf(k1m, 1e-3f);
                    if (lam * lam > 0.98f * k1 * k2) lam = 0.9899f * sqrtf(k1 * k2) * (rho >= 0 ? 1.f : -1.f);
                    s_mu[k] = mu; s_nu[k] = nu; s_k1[k] = k1; s_k2[k] = k2; s_lam[k] = lam;
                    s_lc[k] = logC(tab, nk, mmax, logkmin, dlogk, k1, k2, lam);
                }
                __syncthreads();
            }
            __syncthreads();
        }

        // ---- 计算拟合优度:网格上的 log 似然(用质量加权)
        float ll = 0.f;
        for (int c = tid; c < kCells; c += BS) {
            int ix = c / kBins, iy = c % kBins;
            float cph = s_cph[ix], sph = s_sph[ix], cps = s_cps[iy], sps = s_sps[iy];
            float mx = -3e38f;
            float tmp[16];
            for (int k = 0; k < KK; ++k) {
                float cm = __cosf(s_mu[k]), sm = __sinf(s_mu[k]);
                float cn = __cosf(s_nu[k]), sn = __sinf(s_nu[k]);
                float c1 = cph * cm + sph * sm, q1 = sph * cm - cph * sm;
                float c2 = cps * cn + sps * sn, q2 = sps * cn - cps * sn;
                tmp[k] = __logf(fmaxf(s_w[k], 1e-30f)) + s_lc[k] + s_k1[k] * c1 + s_k2[k] * c2 + s_lam[k] * q1 * q2;
                if (tmp[k] > mx) mx = tmp[k];
            }
            float sum = 0.f;
            for (int k = 0; k < KK; ++k) sum += __expf(tmp[k] - mx);
            float ld = mx + __logf(fmaxf(sum, 1e-30f));
            ll += s_mass[c] * ld;
        }
        ll = block_sum(ll, s_red);
        if (tid == 0) {
            float bic = -2.f * ll * fmaxf(float(recs[ctx].n_obs), 1.f) + float(KK) * 6.f * __logf(fmaxf(float(recs[ctx].n_obs), 2.f));
            if (bic < best_bic) {
                best_bic = bic;
                best_k = KK;
                ParamRec R{};
                R.tier = recs[ctx].tier;
                R.n_obs = recs[ctx].n_obs;
                R.k_eff = uint32_t(KK);
                R.n_samples = recs[ctx].n_samples;
                R.k_post_mean = recs[ctx].k_post_mean;
                R.k_post_sd = recs[ctx].k_post_sd;
                R.loglik_mean = ll;
                R.reserved_f = bic;
                for (int k = 0; k < kKMax; ++k) {
                    float* d = R.comp + k * 8;
                    if (k < KK) {
                        d[0] = s_w[k]; d[1] = s_mu[k]; d[2] = s_nu[k];
                        d[3] = s_k1[k]; d[4] = s_k2[k]; d[5] = s_lam[k]; d[6] = s_lc[k];
                    } else {
                        d[0] = 0.f; d[1] = 0.f; d[2] = 0.f; d[3] = 1.f; d[4] = 1.f; d[5] = 0.f; d[6] = 0.f;
                    }
                    d[7] = 0.f;
                }
                bestrec = R;
            }
        }
        __syncthreads();
    }
    if (tid == 0) recs[ctx] = bestrec;
}

// ================================================================ 采样
__global__ void k_sample(const ParamRec* __restrict__ recs, uint32_t* __restrict__ out,
                         uint64_t m_samples, uint64_t seed, int sweeps) {
    const int ctx = blockIdx.x;
    const int tid = threadIdx.x;
    const int BS = blockDim.x;
    __shared__ uint32_t s_cnt[kCells];

    for (int c = tid; c < kCells; c += BS) s_cnt[c] = 0u;
    __syncthreads();

    ParamRec R = recs[ctx];
    int K = int(R.k_eff);
    if (K < 1) { if (tid == 0) R.tier = R.tier; return; }
    if (K > kKMax) K = kKMax;

    // 归一化权重
    float cw[kKMax];
    float wsum = 0.f;
    for (int k = 0; k < K; ++k) { cw[k] = fmaxf(R.comp[k * 8 + 0], 0.f); wsum += cw[k]; }
    if (wsum <= 0.f) { cw[0] = 1.f; wsum = 1.f; }
    for (int k = 0; k < K; ++k) cw[k] /= wsum;

    uint64_t per = m_samples / uint64_t(BS) + 1;
    for (uint64_t s = 0; s < per; ++s) {
        uint64_t idx = uint64_t(tid) * per + s;
        if (idx >= m_samples) continue;
        Rng r = make_rng(seed, ctx, uint32_t(tid), uint32_t(s), 0x5EEDu);
        float u = u01(r);
        int pick = K - 1;
        float acc = 0.f;
        for (int k = 0; k < K; ++k) { acc += cw[k]; if (u <= acc) { pick = k; break; } }
        float mu = R.comp[pick * 8 + 1], nu = R.comp[pick * 8 + 2];
        float k1 = R.comp[pick * 8 + 3], k2 = R.comp[pick * 8 + 4], lam = R.comp[pick * 8 + 5];
        float ph = mu, ps = nu;
        for (int it = 0; it < sweeps; ++it) {
            float s2 = __sinf(ps - nu);
            float b1 = lam * s2;
            ph = sample_vm(mu + atan2f(b1, k1), sqrtf(k1 * k1 + b1 * b1), r);
            float s1 = __sinf(ph - mu);
            float b2 = lam * s1;
            ps = sample_vm(nu + atan2f(b2, k2), sqrtf(k2 * k2 + b2 * b2), r);
        }
        ph = fmodf(ph + 3.14159265358979f, 6.28318530718f);
        if (ph < 0.f) ph += 6.28318530718f;
        ps = fmodf(ps + 3.14159265358979f, 6.28318530718f);
        if (ps < 0.f) ps += 6.28318530718f;
        int ix = int(ph / (6.28318530718f / float(kBins)));
        int iy = int(ps / (6.28318530718f / float(kBins)));
        if (ix < 0) ix = 0; if (ix >= kBins) ix = kBins - 1;
        if (iy < 0) iy = 0; if (iy >= kBins) iy = kBins - 1;
        atomicAdd(&s_cnt[ix * kBins + iy], 1u);
    }
    __syncthreads();
    for (int c = tid; c < kCells; c += BS) out[size_t(ctx) * kCells + c] = s_cnt[c];
}

}  // namespace

// ================================================================ 宿主
bool gpu_info(GpuInfo& out, std::string& err) {
    int dev = 0;
    cudaError_t e = cudaGetDevice(&dev);
    if (e != cudaSuccess) { err = std::string("cudaGetDevice: ") + estr(e); return false; }
    cudaDeviceProp p{};
    e = cudaGetDeviceProperties(&p, dev);
    if (e != cudaSuccess) { err = std::string("props: ") + estr(e); return false; }
    std::snprintf(out.name, sizeof(out.name), "%s", p.name);
    out.vram = p.totalGlobalMem;
    out.cc_major = p.major; out.cc_minor = p.minor;
    out.sm_count = p.multiProcessorCount;
    return true;
}

bool gpu_init(const BesselTable& tab, std::string& err) {
    build_lbinom();
    g.nk = tab.nk; g.mmax = tab.mmax;
    g.logkmin = float(std::log(tab.kmin));
    g.dlogk = float(tab.dlogk);
    cudaError_t e;
    e = cudaMalloc(&g.tab, tab.v.size() * sizeof(float));
    if (e != cudaSuccess) { err = "cudaMalloc bessel"; return false; }
    e = cudaMemcpy(g.tab, tab.v.data(), tab.v.size() * sizeof(float), cudaMemcpyHostToDevice);
    if (e != cudaSuccess) { err = "memcpy bessel"; return false; }

    e = cudaMalloc(&g.d_hist, size_t(kNctx) * kCells * sizeof(uint32_t));
    if (e != cudaSuccess) { err = "cudaMalloc hist"; return false; }
    e = cudaMalloc(&g.d_recs, sizeof(ParamRec) * size_t(kNctx));
    if (e != cudaSuccess) { err = "cudaMalloc recs"; return false; }
    e = cudaMalloc(&g.d_grid, size_t(kNctx) * kCells * sizeof(float));
    if (e != cudaSuccess) { err = "cudaMalloc grid"; return false; }
    g.ready = true;
    return true;
}

void gpu_shutdown() {
    cudaFree(g.tab); cudaFree(g.d_hist); cudaFree(g.d_recs); cudaFree(g.d_grid);
    cudaFree(g.d_start); cudaFree(g.d_order);
    cudaFree(g.d_cphi); cudaFree(g.d_sphi); cudaFree(g.d_cpsi); cudaFree(g.d_spsi);
    cudaFree(g.d_bx); cudaFree(g.d_by); cudaFree(g.d_label);
    g = Dev{};
}

static bool upload_points(const PointSet& ps, std::string& err) {
    size_t n = ps.size();
    if (n > g.cap_points) {
        cudaFree(g.d_cphi); cudaFree(g.d_sphi); cudaFree(g.d_cpsi); cudaFree(g.d_spsi);
        cudaFree(g.d_bx); cudaFree(g.d_by); cudaFree(g.d_label);
        cudaFree(g.d_start); cudaFree(g.d_order);
        g.d_cphi = g.d_sphi = g.d_cpsi = g.d_spsi = nullptr;
        g.d_bx = g.d_by = nullptr; g.d_label = nullptr;
        g.d_start = g.d_order = nullptr;
        g.cap_points = 0;
        if (cudaMalloc(&g.d_cphi, n * sizeof(float)) != cudaSuccess) { err = "malloc cphi"; return false; }
        if (cudaMalloc(&g.d_sphi, n * sizeof(float)) != cudaSuccess) { err = "malloc sphi"; return false; }
        if (cudaMalloc(&g.d_cpsi, n * sizeof(float)) != cudaSuccess) { err = "malloc cpsi"; return false; }
        if (cudaMalloc(&g.d_spsi, n * sizeof(float)) != cudaSuccess) { err = "malloc spsi"; return false; }
        if (cudaMalloc(&g.d_bx, n) != cudaSuccess) { err = "malloc bx"; return false; }
        if (cudaMalloc(&g.d_by, n) != cudaSuccess) { err = "malloc by"; return false; }
        if (cudaMalloc(&g.d_label, n * sizeof(int)) != cudaSuccess) { err = "malloc label"; return false; }
        if (cudaMalloc(&g.d_start, size_t(kNctx + 1) * sizeof(uint32_t)) != cudaSuccess) { err = "malloc start"; return false; }
        if (cudaMalloc(&g.d_order, n * sizeof(uint32_t)) != cudaSuccess) { err = "malloc order"; return false; }
        g.cap_points = n;
    }
    #define UPL(dst, src, cnt, esz) \
        if (cudaMemcpy(dst, src, size_t(cnt) * esz, cudaMemcpyHostToDevice) != cudaSuccess) { err = "memcpy " #dst; return false; }
    UPL(g.d_cphi, ps.cphi.data(), n, sizeof(float))
    UPL(g.d_sphi, ps.sphi.data(), n, sizeof(float))
    UPL(g.d_cpsi, ps.cpsi.data(), n, sizeof(float))
    UPL(g.d_spsi, ps.spsi.data(), n, sizeof(float))
    UPL(g.d_bx, ps.bx.data(), n, 1)
    UPL(g.d_by, ps.by.data(), n, 1)
    UPL(g.d_start, ps.start.data(), kNctx + 1, sizeof(uint32_t))
    UPL(g.d_order, ps.order.data(), n, sizeof(uint32_t))
    #undef UPL
    // 重排 cos/sin 到 CSR 顺序,便于按簇连续访问
    return true;
}

bool gpu_histogram(const PointSet& ps, std::vector<uint32_t>& counts, std::string& err) {
    if (cudaMemset(g.d_hist, 0, size_t(kNctx) * kCells * sizeof(uint32_t)) != cudaSuccess) {
        err = "memset hist"; return false;
    }
    int bs = 256;
    int nb = int((ps.n_points + bs - 1) / bs);
    if (nb > 65535) nb = 65535;
    k_hist<<<nb, bs>>>(g.d_bx, g.d_by, ps.ctx.data(), ps.n_points, g.d_hist);
    if (cudaDeviceSynchronize() != cudaSuccess) { err = "k_hist"; return false; }
    counts.resize(size_t(kNctx) * kCells);
    if (cudaMemcpy(counts.data(), g.d_hist, counts.size() * sizeof(uint32_t),
                   cudaMemcpyDeviceToHost) != cudaSuccess) { err = "copy hist"; return false; }
    return true;
}

bool gpu_fit(const Config& cfg, const PointSet& ps, const CtxStat& st,
             std::vector<ParamRec>& recs, std::vector<float>& pred_grid, std::string& err) {
    if (!upload_points(ps, err)) return false;

    // 分级 + 按代价降序排布(长任务先发)
    struct Item { uint32_t ctx; uint8_t tier; uint32_t n; };
    std::vector<Item> items;
    items.reserve(kNctx);
    uint64_t nz = 0;
    for (int c = 0; c < kNctx; ++c) {
        uint32_t n = st.n[size_t(c)];
        uint8_t tier = 0;
        if (n == 0) tier = 0;
        else if (int(n) >= cfg.tier_a_n) tier = 3;
        else if (int(n) >= cfg.tier_b_n) tier = 2;
        else tier = 1;
        if (n > 0) { items.push_back({uint32_t(c), tier, n}); ++nz; }
    }
    std::sort(items.begin(), items.end(), [](const Item& a, const Item& b) {
        if (a.tier != b.tier) return a.tier > b.tier;
        return a.n > b.n;
    });

    std::vector<uint32_t> ctx_list(items.size());
    std::vector<uint8_t>  tier_list(items.size());
    for (size_t i = 0; i < items.size(); ++i) { ctx_list[i] = items[i].ctx; tier_list[i] = items[i].tier; }

    uint32_t *d_ctx = nullptr, *d_n = nullptr;
    uint8_t* d_tier = nullptr;
    if (cudaMalloc(&d_ctx, ctx_list.size() * sizeof(uint32_t)) != cudaSuccess) { err = "malloc ctxlist"; return false; }
    if (cudaMalloc(&d_tier, tier_list.size()) != cudaSuccess) { err = "malloc tierlist"; return false; }
    cudaMemcpy(d_ctx, ctx_list.data(), ctx_list.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_tier, tier_list.data(), tier_list.size(), cudaMemcpyHostToDevice);
    (void)d_n;

    // 清零网格 + 初始化记录
    cudaMemset(g.d_grid, 0, size_t(kNctx) * kCells * sizeof(float));
    recs.assign(size_t(kNctx), ParamRec{});
    for (int c = 0; c < kNctx; ++c) {
        recs[size_t(c)].tier = 0;
        recs[size_t(c)].n_obs = st.n[size_t(c)];
        recs[size_t(c)].k_eff = 0;
        recs[size_t(c)].n_samples = 0;
        for (int k = 0; k < kKMax; ++k) {
            float* d = recs[size_t(c)].comp + k * 8;
            d[0] = 0.f; d[1] = 0.f; d[2] = 0.f; d[3] = 1.f; d[4] = 1.f; d[5] = 0.f; d[6] = 0.f; d[7] = 0.f;
        }
        if (st.n[size_t(c)] == 0) {
            // D 级:先验预测(近似均匀)
            float* d = recs[size_t(c)].comp;
            d[0] = 1.f; d[1] = 0.f; d[2] = 0.f;
            d[3] = cfg.k10; d[4] = cfg.k20; d[5] = cfg.lam0;
            d[6] = float(logC_exact(cfg.k10, cfg.k20, cfg.lam0));
            recs[size_t(c)].k_eff = 1;
        }
    }
    cudaMemcpy(g.d_recs, recs.data(), sizeof(ParamRec) * size_t(kNctx), cudaMemcpyHostToDevice);

    FitArgs A{};
    A.start = g.d_start; A.order = g.d_order;
    A.cphi = g.d_cphi; A.sphi = g.d_sphi; A.cpsi = g.d_cpsi; A.spsi = g.d_spsi;
    A.label = g.d_label;
    A.tab = g.tab; A.nk = g.nk; A.mmax = g.mmax; A.logkmin = g.logkmin; A.dlogk = g.dlogk;
    A.ctx_list = d_ctx; A.tier_list = d_tier;
    A.tau0 = cfg.tau0; A.k10 = cfg.k10; A.k20 = cfg.k20; A.lam0 = cfg.lam0;
    A.nu0 = cfg.wish_nu0; A.b0 = cfg.wish_b0; A.k_aux = cfg.k_aux;
    A.seed = cfg.seed;
    A.iters_a = cfg.a_iters; A.burn_a = cfg.a_burn; A.thin_a = cfg.a_thin; A.chains_a = cfg.a_chains;
    A.iters_b = cfg.b_iters; A.burn_b = cfg.b_burn; A.thin_b = cfg.b_thin; A.chains_b = cfg.b_chains;
    A.iters_c = cfg.c_iters; A.burn_c = cfg.c_burn; A.thin_c = cfg.c_thin; A.chains_c = cfg.c_chains;
    A.recs = g.d_recs; A.grid = g.d_grid;

    if (!items.empty()) {
        int bs = 128;
        k_mcmc<<<int(items.size()), bs>>>(A);
        if (cudaDeviceSynchronize() != cudaSuccess) { err = "k_mcmc"; return false; }
    }
    // D 级网格:先验预测(在宿主填)
    {
        std::vector<float> gbuf(size_t(kNctx) * kCells, 0.f);
        cudaMemcpy(gbuf.data(), g.d_grid, gbuf.size() * sizeof(float), cudaMemcpyDeviceToHost);
        double lc = logC_exact(cfg.k10, cfg.k20, cfg.lam0);
        std::vector<double> lcb(kCells);
        for (int c = 0; c < kCells; ++c) {
            int ix = c / kBins, iy = c % kBins;
            double ph = bin_center(ix), ps = bin_center(iy);
            lcb[c] = lc + cfg.k10 * std::cos(ph) + cfg.k20 * std::cos(ps)
                   + cfg.lam0 * std::sin(ph) * std::sin(ps);
        }
        for (int c = 0; c < kNctx; ++c) {
            if (st.n[size_t(c)] != 0) continue;
            for (int j = 0; j < kCells; ++j) gbuf[size_t(c) * kCells + j] = float(std::exp(lcb[j]));
        }
        cudaMemcpy(g.d_grid, gbuf.data(), gbuf.size() * sizeof(float), cudaMemcpyHostToDevice);
    }

    // 压缩
    k_compress<<<kNctx, 128>>>(g.d_grid, g.d_recs, g.tab, g.nk, g.mmax, g.logkmin, g.dlogk,
                               cfg.k_max_comp, 40);
    if (cudaDeviceSynchronize() != cudaSuccess) { err = "k_compress"; return false; }

    cudaMemcpy(recs.data(), g.d_recs, sizeof(ParamRec) * size_t(kNctx), cudaMemcpyDeviceToHost);
    pred_grid.resize(size_t(kNctx) * kCells);
    cudaMemcpy(pred_grid.data(), g.d_grid, pred_grid.size() * sizeof(float), cudaMemcpyDeviceToHost);

    // 归一化预测网格(密度 -> 格概率)
    for (int c = 0; c < kNctx; ++c) {
        float* s = pred_grid.data() + size_t(c) * kCells;
        double tot = 0.0;
        for (int i = 0; i < kCells; ++i) { s[i] = float(double(s[i]) * kCellArea); tot += s[i]; }
        float inv = tot > 0 ? float(1.0 / tot) : 0.f;
        for (int i = 0; i < kCells; ++i) s[i] *= inv;
    }

    cudaFree(d_ctx); cudaFree(d_tier);
    return true;
}

bool gpu_sample(const Config& cfg, const std::vector<ParamRec>& recs,
                std::vector<uint32_t>& counts, std::string& err) {
    cudaMemcpy(g.d_recs, recs.data(), sizeof(ParamRec) * size_t(kNctx), cudaMemcpyHostToDevice);
    cudaMemset(g.d_hist, 0, size_t(kNctx) * kCells * sizeof(uint32_t));
    k_sample<<<kNctx, 128>>>(g.d_recs, g.d_hist, cfg.m_samples, cfg.seed ^ 0xABCDEFull, 6);
    if (cudaDeviceSynchronize() != cudaSuccess) { err = "k_sample"; return false; }
    counts.resize(size_t(kNctx) * kCells);
    cudaMemcpy(counts.data(), g.d_hist, counts.size() * sizeof(uint32_t), cudaMemcpyDeviceToHost);
    return true;
}

double host_logC(const BesselTable&, double k1, double k2, double lam) {
    return logC_exact(k1, k2, lam);
}

}  // namespace vd



