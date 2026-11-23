// von_dpmm / main.cpp — 编排:读 bin -> CSR -> 直方图(staic_box) -> MCMC(von) -> 采样(von_box)
#include "model.h"
#include "gpu.h"

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <filesystem>
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>

namespace fs = std::filesystem;
using namespace vd;

static double now_s() {
    using clk = std::chrono::steady_clock;
    return std::chrono::duration<double>(clk::now().time_since_epoch()).count();
}

static void say(const std::string& s) { std::printf("%s\n", s.c_str()); std::fflush(stdout); }

static void progress(const char* tag, double frac, const std::string& extra) {
    std::printf("\r[%5.1f%%] %-10s %s", frac * 100.0, tag, extra.c_str());
    std::fflush(stdout);
}

int main(int argc, char** argv) {
    SetConsoleOutputCP(65001);
    std::string cfg_path = "config.json";
    bool selftest = false, verify = false, stage_hist = true, stage_fit = true, stage_sample = true;
    bool do_recompress = false;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if ((a == "--config" || a == "-c") && i + 1 < argc) cfg_path = argv[++i];
        else if (a == "--selftest") selftest = true;
        else if (a == "--verify") verify = true;
        else if (a == "--no-hist") stage_hist = false;
        else if (a == "--no-fit") stage_fit = false;
        else if (a == "--no-sample") stage_sample = false;
        else if (a == "--only-hist") { stage_fit = false; stage_sample = false; }
        else if (a == "--recompress") { do_recompress = true; stage_hist = false; stage_fit = false; stage_sample = false; }
        else if (a == "--help" || a == "-h") {
            say("von_dpmm [--config path] [--selftest] [--verify] [--only-hist] [--no-fit] [--no-sample]");
            return 0;
        }
    }

    Config cfg;
    std::string err;
    if (!load_config(cfg_path, cfg, err)) { say("配置错误: " + err); return 2; }

    // ------------------------------------------------------------ 后端选择
    GpuInfo gi;
    bool use_gpu = false;
    {
        std::string gerr;
        if (gpu_info(gi, gerr)) {
            say(fmt("设备: %s  sm_%d%d  SM=%d  VRAM=%.0f MB", gi.name, gi.cc_major, gi.cc_minor,
                    gi.sm_count, double(gi.vram) / 1048576.0));
            use_gpu = true;
        } else {
            say("GPU 不可用: " + gerr);
            say("  -> 退回 CPU 后端(同算法,多线程)");
        }
    }
    say(fmt("配置: BIN=%s", cfg.bin_dir.c_str()));
    say(fmt("      OUT=%s", cfg.out_root.c_str()));
    say(fmt("      TIER A>=%d  B>=%d  |  A %d/%d/%d/%d链  B %d/%d/%d/%d链  C %d/%d/%d/%d链",
            cfg.tier_a_n, cfg.tier_b_n,
            cfg.a_iters, cfg.a_burn, cfg.a_thin, cfg.a_chains,
            cfg.b_iters, cfg.b_burn, cfg.b_thin, cfg.b_chains,
            cfg.c_iters, cfg.c_burn, cfg.c_thin, cfg.c_chains));
    say(fmt("      tau0=%.3f  H1(k=%.4f,%.4f)  H2(W_%d, B0=%.3f)  K_AUX=%d  SEED=%llu",
            double(cfg.tau0), double(cfg.k10), double(cfg.k20), cfg.wish_nu0,
            double(cfg.wish_b0), cfg.k_aux, (unsigned long long)cfg.seed));

    // ------------------------------------------------------------ Bessel 表
    double t0 = now_s();
    BesselTable tab;
    build_bessel_table(tab);
    say(fmt("Bessel 表: %d x %d  (log k in [%.3g,%.3g])  建表 %.2f s",
            tab.nk, tab.mmax + 1, tab.kmin, tab.kmax, now_s() - t0));

    if (selftest) {
        // 自检 1:Bessel 表 vs 直算
        double mx = 0.0;
        for (int m = 0; m <= 8; ++m) {
            for (double k : {0.01, 0.1, 1.0, 5.0, 20.0, 50.0, 200.0, 600.0}) {
                double exact = logIe_exact(m, k);
                double lk = std::log(k);
                double t = (lk - std::log(tab.kmin)) / tab.dlogk;
                int i0 = int(t); double f = t - i0;
                if (i0 < 0) { i0 = 0; f = 0; }
                if (i0 >= tab.nk - 1) { i0 = tab.nk - 2; f = 1; }
                const float* r0 = tab.v.data() + size_t(i0) * (tab.mmax + 1) + m;
                const float* r1 = tab.v.data() + size_t(i0 + 1) * (tab.mmax + 1) + m;
                double interp = r0[0] * (1.0 - f) + r1[0] * f;
                if (exact > -1e29) mx = std::max(mx, std::fabs(exact - interp));
            }
        }
        say(fmt("自检 logIe 插值最大偏差: %.3e", mx));

        // 自检 2:lambda=0 时 C 应等于 1/(4pi^2 I0(k1) I0(k2))
        double e2 = 0.0;
        for (double k1 : {1.0, 10.0, 100.0}) {
            for (double k2 : {1.0, 10.0, 100.0}) {
                double c1 = logC_exact(k1, k2, 0.0);
                double ref = -std::log(4.0 * kPi * kPi) - std::log(std::exp(logIe_exact(0, k1)) * std::exp(k1))
                           - std::log(std::exp(logIe_exact(0, k2)) * std::exp(k2));
                e2 = std::max(e2, std::fabs(c1 - ref));
            }
        }
        say(fmt("自检 logC(lam=0) 与 1/(4pi^2 I0 I0) 偏差: %.3e", e2));

        // 自检 3:归一化数值积分
        {
            double k1 = 20, k2 = 30, lam = 8;
            double lc = logC_exact(k1, k2, lam);
            int NG = 720;
            double d = kTwoPi / NG;
            double sum = 0.0;
            for (int i = 0; i < NG; ++i) {
                double ph = -kPi + (i + 0.5) * d;
                for (int j = 0; j < NG; ++j) {
                    double ps = -kPi + (j + 0.5) * d;
                    sum += std::exp(lc + k1 * std::cos(ph) + k2 * std::cos(ps)
                                    + lam * std::sin(ph) * std::sin(ps)) * d * d;
                }
            }
            say(fmt("自检 归一化积分 (k1=20,k2=30,lam=8): %.6f  (应=1)", sum));
        }
        say("自检完成");
        return 0;
    }

    if (use_gpu && !gpu_init(tab, err)) { say("GPU 初始化失败: " + err + " -> CPU"); use_gpu = false; }
say(fmt("后端: %s", use_gpu ? "CUDA" : "CPU"));

    // ------------------------------------------------------------ 读数据
    t0 = now_s();
    PointSet ps;
    if (!load_points(cfg, ps, err)) { say("读数据失败: " + err); return 4; }
    say(fmt("扫描: 段 %llu, 残基 %llu, 可用点 %llu  (边界丢 %llu / NaN %llu / 范数 %llu / 非标准 %llu)  %.2f s",
            (unsigned long long)ps.n_seg, (unsigned long long)ps.n_res,
            (unsigned long long)ps.n_points, (unsigned long long)ps.drop_edge,
            (unsigned long long)ps.drop_nan, (unsigned long long)ps.drop_norm,
            (unsigned long long)ps.drop_kine, now_s() - t0));

    t0 = now_s();
    CtxStat st;
    if (!build_csr(cfg, ps, err, st)) { say("CSR 失败: " + err); return 4; }
    say(fmt("CSR 建好  %.2f s", now_s() - t0));

    // ------------------------------------------------------------ 分级统计
    int nA = 0, nB = 0, nC = 0, nD = 0;
    for (int c = 0; c < kNctx; ++c) {
        int n = int(st.n[size_t(c)]);
        if (n == 0) ++nD;
        else if (n >= cfg.tier_a_n) ++nA;
        else if (n >= cfg.tier_b_n) ++nB;
        else ++nC;
    }
    say(fmt("分级: A=%d  B=%d  C=%d  D(空)=%d   合计 %d", nA, nB, nC, nD, kNctx));
    std::vector<std::pair<int, uint32_t>> top;
    for (int c = 0; c < kNctx; ++c) top.push_back({int(st.n[size_t(c)]), uint32_t(c)});
    std::sort(top.rbegin(), top.rend());
    say("样本最多的 12 个上下文 (m-1,m0,m+1 : N):");
    for (int i = 0; i < 12 && i < int(top.size()); ++i) {
        uint32_t c = top[i].second;
        int a = int(c) / (kNaa * kNaa), b = (int(c) / kNaa) % kNaa, cc = int(c) % kNaa;
        say(fmt("   %c%c%c : %d", kAaOrder[a], kAaOrder[b], kAaOrder[cc], top[i].first));
    }
    {
        std::vector<int> ns;
        for (int c = 0; c < kNctx; ++c) if (st.n[size_t(c)] > 0) ns.push_back(int(st.n[size_t(c)]));
        std::sort(ns.begin(), ns.end());
        auto q = [&](double p) { return ns.empty() ? 0 : ns[size_t(p * (ns.size() - 1))]; };
        say(fmt("N 分位: min=%d  p10=%d  p25=%d  中位=%d  p75=%d  p90=%d  max=%d  (非空 %zu 个)",
                ns.empty() ? 0 : ns.front(), q(0.10), q(0.25), q(0.50), q(0.75), q(0.90),
                ns.empty() ? 0 : ns.back(), ns.size()));
    }

    fs::create_directories(cfg.out_root + "\\von");
    fs::create_directories(cfg.out_root + "\\von_box");
    fs::create_directories(cfg.out_root + "\\staic_box");
    fs::create_directories(cfg.out_root + "\\img");

    // ------------------------------------------------------------ 重压缩(不重跑 MCMC)
    if (do_recompress) {
        t0 = now_s();
        std::string pp = cfg.out_root + "\\von\\von_params.bin";
        std::string pg = cfg.out_root + "\\von\\von_pred_grid.bin";
        FILE* fp = fopen(pp.c_str(), "rb");
        if (!fp) { say("缺少 " + pp); return 9; }
        ParamHeader ph{};
        if (fread(&ph, sizeof(ph), 1, fp) != 1 || std::memcmp(ph.magic, "VPM1", 4) != 0) {
            say("参数文件头不对"); fclose(fp); return 9;
        }
        std::vector<ParamRec> rr;
        rr.resize(size_t(kNctx));
        if (fread(rr.data(), sizeof(ParamRec), size_t(kNctx), fp) != size_t(kNctx)) {
            say("参数文件被截断"); fclose(fp); return 9;
        }
        fclose(fp);

        std::vector<float> gg;
        gg.resize(size_t(kNctx) * kCells);
        fp = fopen(pg.c_str(), "rb");
        if (!fp) { say("缺少 " + pg); return 9; }
        if (fseek(fp, 40, SEEK_SET) != 0) { fclose(fp); return 9; }
        if (fread(gg.data(), sizeof(float), gg.size(), fp) != gg.size()) {
            say("预测网格被截断"); fclose(fp); return 9;
        }
        fclose(fp);

        double mtv = 0.0;
        say("按式(7) 网格重做压缩(K' 由 BIC 选,初值=质量分位,无 MCMC 种子)…");
        if (!cpu_recompress(cfg, rr, gg.data(), &mtv, err)) { say("重压缩失败: " + err); return 9; }
        say(fmt("重压缩完成  平均 TV(压缩模型网格, 式7网格) = %.4f   %.1f s", mtv, now_s() - t0));
        int ksum = 0;
        for (int c = 0; c < kNctx; ++c) ksum += int(rr[size_t(c)].k_eff);
        say(fmt("  平均 K' = %.2f", double(ksum) / kNctx));
        if (!write_params(pp, rr.data(), ps.n_points, err)) { say("写参数失败: " + err); return 9; }
        say("  已写回 " + pp);
        if (use_gpu) gpu_shutdown();
        say("完成");
        return 0;
    }

    // ------------------------------------------------------------ staic_box
    std::vector<uint32_t> hist;
    if (stage_hist) {
        t0 = now_s();
        bool okh = use_gpu ? gpu_histogram(ps, hist, err) : cpu_histogram(ps, hist);
        if (!okh) { say("直方图失败: " + err); return 5; }
        uint64_t tot = 0;
        for (uint32_t v : hist) tot += v;
        say(fmt("staic_box 直方图: 总计数 %llu (点 %llu)  差异 %lld  %.2f s",
                (unsigned long long)tot, (unsigned long long)ps.n_points,
                (long long)tot - (long long)ps.n_points, now_s() - t0));
        std::string p = cfg.out_root + "\\staic_box\\staic_box.bin";
        if (!write_staic_box(p, hist.data(), ps.n_points, err)) { say("写 staic_box 失败: " + err); return 6; }
        say(fmt("  写出 %s  (%.1f MB)", p.c_str(),
                double(size_t(kNctx) * kCells * 8 + sizeof(BoxHeader)) / 1048576.0));
    }

    // ------------------------------------------------------------ von (MCMC)
    std::vector<ParamRec> recs;
    std::vector<float> pred;
    if (stage_fit) {
        t0 = now_s();
        say("DPMM-MCMC 开始(单次核启动,无中间进度)…");
        bool okf = use_gpu ? gpu_fit(cfg, ps, st, recs, pred, err)
                           : cpu_fit(cfg, ps, st, recs, pred, err);
        if (!okf) { say("MCMC 失败: " + err); return 7; }
        say(fmt("MCMC + 压缩完成  %.1f s", now_s() - t0));
        std::string pp = cfg.out_root + "\\von\\von_params.bin";
        if (!write_params(pp, recs.data(), ps.n_points, err)) { say("写参数失败: " + err); return 6; }
        std::string pg = cfg.out_root + "\\von\\von_pred_grid.bin";
        if (!write_pred_grid(pg, pred.data(), err)) { say("写预测网格失败: " + err); return 6; }
        say(fmt("  写出 %s (%.2f MB) 与 %s (%.1f MB)", pp.c_str(),
                double(size_t(kNctx) * sizeof(ParamRec)) / 1048576.0, pg.c_str(),
                double(size_t(kNctx) * kCells * 4) / 1048576.0));
        // 参数统计
        int ksum = 0, kmaxseen = 0;
        for (int c = 0; c < kNctx; ++c) { ksum += int(recs[size_t(c)].k_eff); kmaxseen = std::max(kmaxseen, int(recs[size_t(c)].k_eff)); }
        say(fmt("  压缩成分: 平均 K'=%.2f  最大 %d", double(ksum) / kNctx, kmaxseen));
        double sd = 0.0; int cnt = 0;
        for (int c = 0; c < kNctx; ++c) if (recs[size_t(c)].n_obs > 0) { sd += recs[size_t(c)].k_post_mean; ++cnt; }
        if (cnt) say(fmt("  MCMC 后验占用簇数均值: %.2f (非空上下文 %d 个)", sd / cnt, cnt));
    }

    // ------------------------------------------------------------ 只重采样:从磁盘读参数
    if (stage_sample && !stage_fit && recs.empty()) {
        std::string pp = cfg.out_root + "\\von\\von_params.bin";
        FILE* fp = fopen(pp.c_str(), "rb");
        if (!fp) { say("缺少 " + pp); return 9; }
        ParamHeader ph{};
        if (fread(&ph, sizeof(ph), 1, fp) != 1 || std::memcmp(ph.magic, "VPM1", 4) != 0) {
            say("参数文件头不对"); fclose(fp); return 9;
        }
        recs.resize(size_t(kNctx));
        if (fread(recs.data(), sizeof(ParamRec), size_t(kNctx), fp) != size_t(kNctx)) {
            say("参数文件被截断"); fclose(fp); return 9;
        }
        fclose(fp);
        say("已从磁盘读入 von_params.bin(跳过 MCMC)");
    }

    // ------------------------------------------------------------ von_box
    if (stage_sample && !recs.empty()) {
        t0 = now_s();
        std::vector<uint32_t> sc;
        bool oks = use_gpu ? gpu_sample(cfg, recs, sc, err) : cpu_sample(cfg, recs, sc, err);
        if (!oks) { say("采样失败: " + err); return 8; }
        uint64_t tot = 0;
        for (uint32_t v : sc) tot += v;
        say(fmt("von_box 采样: 总样本 %llu (期望 %llu)  %.2f s",
                (unsigned long long)tot, (unsigned long long)(uint64_t(kNctx) * cfg.m_samples),
                now_s() - t0));
        std::string p = cfg.out_root + "\\von_box\\von_box.bin";
        if (!write_von_box(p, sc.data(), tot, err)) { say("写 von_box 失败: " + err); return 6; }
        say(fmt("  写出 %s", p.c_str()));

        // 一致性:von_box vs 预测网格 vs 经验
        if (!hist.empty() && !pred.empty()) {
            double tv_pred = 0.0, tv_emp = 0.0;
            int cnt = 0;
            double worst = -1; int worst_ctx = -1;
            for (int c = 0; c < kNctx; ++c) {
                if (st.n[size_t(c)] < 200) continue;
                double t1 = 0.0, t2 = 0.0;
                double s_sc = 0, s_pr = 0, s_hi = 0;
                for (int i = 0; i < kCells; ++i) {
                    s_sc += sc[size_t(c) * kCells + i];
                    s_pr += pred[size_t(c) * kCells + i];
                    s_hi += hist[size_t(c) * kCells + i];
                }
                if (s_sc <= 0 || s_pr <= 0 || s_hi <= 0) continue;
                for (int i = 0; i < kCells; ++i) {
                    double a = sc[size_t(c) * kCells + i] / s_sc;
                    double b = pred[size_t(c) * kCells + i] / s_pr;
                    double e = hist[size_t(c) * kCells + i] / s_hi;
                    t1 += std::fabs(a - b);
                    t2 += std::fabs(b - e);
                }
                tv_pred += t1; tv_emp += t2; ++cnt;
                if (t2 > worst) { worst = t2; worst_ctx = c; }
            }
            if (cnt) {
                say(fmt("TV(von_box vs 式7网格) 均值 %.4f   TV(式7网格 vs 经验) 均值 %.4f  (N>=200 的 %d 个)",
                        tv_pred / cnt, tv_emp / cnt, cnt));
                say(fmt("最差上下文 %d 的 TV(vs 经验) = %.4f", worst_ctx, worst));
            }
        }
    }

    if (verify) {
        say("verify: 检查输出文件头");
        for (const char* f : { "\\staic_box\\staic_box.bin", "\\von\\von_params.bin",
                               "\\von\\von_pred_grid.bin", "\\von_box\\von_box.bin" }) {
            std::string p = cfg.out_root + f;
            FILE* fp = fopen(p.c_str(), "rb");
            if (!fp) { say(fmt("  缺失 %s", p.c_str())); continue; }
            char m[5] = {0};
            fread(m, 1, 4, fp);
            fseek(fp, 0, SEEK_END);
            long sz = ftell(fp);
            fclose(fp);
            say(fmt("  %-34s magic=%-4s  %ld bytes", f, m, sz));
        }
    }

    if (use_gpu) gpu_shutdown();
    say("完成");
    return 0;
}






