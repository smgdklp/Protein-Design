// von_dpmm / config.cpp — 容错 JSON 配置(键名大小写不敏感,漏逗号也能读)
#include "model.h"

#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include <cstring>
#include <cctype>
#include <fstream>
#include <sstream>
#include <unordered_map>

namespace vd {

static std::string lower(std::string s) {
    for (char& c : s) c = char(std::tolower((unsigned char)c));
    return s;
}

static std::string trim(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

std::string fmt(const char* f, ...) {
    char buf[4096];
    va_list ap;
    va_start(ap, f);
    vsnprintf(buf, sizeof(buf), f, ap);
    va_end(ap);
    return std::string(buf);
}

bool load_config(const std::string& path, Config& cfg, std::string& err) {
    std::ifstream in(path, std::ios::binary);
    if (!in) { err = "cannot open config: " + path; return false; }
    std::stringstream ss;
    ss << in.rdbuf();
    std::string text = ss.str();

    std::unordered_map<std::string, std::string> kv;
    size_t i = 0;
    while (i < text.size()) {
        size_t q1 = text.find('"', i);
        if (q1 == std::string::npos) break;
        size_t q2 = text.find('"', q1 + 1);
        if (q2 == std::string::npos) break;
        std::string key = lower(trim(text.substr(q1 + 1, q2 - q1 - 1)));
        size_t colon = text.find(':', q2 + 1);
        if (colon == std::string::npos) break;
        size_t p = colon + 1;
        while (p < text.size() && std::isspace((unsigned char)text[p])) ++p;
        std::string val;
        if (p < text.size() && text[p] == '"') {
            size_t e = p + 1;
            while (e < text.size() && text[e] != '"') ++e;
            val = text.substr(p + 1, e - p - 1);
            i = e + 1;
        } else if (p < text.size() && (text[p] == 't' || text[p] == 'f')) {
            size_t e = p;
            while (e < text.size() && std::isalpha((unsigned char)text[e])) ++e;
            val = text.substr(p, e - p);
            i = e;
        } else {
            size_t e = p;
            while (e < text.size() && text[e] != ',' && text[e] != '}' && text[e] != '\n') ++e;
            val = trim(text.substr(p, e - p));
            i = e;
        }
        if (!key.empty()) kv[key] = val;
    }

    auto gets = [&](std::initializer_list<const char*> names, std::string& out) -> bool {
        for (const char* n : names) {
            auto it = kv.find(lower(n));
            if (it != kv.end()) { out = it->second; return true; }
        }
        return false;
    };
    auto geti = [&](std::initializer_list<const char*> names, int& out) {
        std::string s;
        if (gets(names, s) && !s.empty()) out = std::atoi(s.c_str());
    };
    auto getu64 = [&](std::initializer_list<const char*> names, uint64_t& out) {
        std::string s;
        if (gets(names, s) && !s.empty()) out = std::strtoull(s.c_str(), nullptr, 10);
    };
    auto getf = [&](std::initializer_list<const char*> names, float& out) {
        std::string s;
        if (gets(names, s) && !s.empty()) out = float(std::atof(s.c_str()));
    };
    auto getb = [&](std::initializer_list<const char*> names, bool& out) {
        std::string s;
        if (gets(names, s)) {
            std::string t = lower(s);
            out = (t == "true" || t == "1" || t == "yes");
        }
    };

    gets({"BIN_DIR", "BIN", "PBD_BIN"}, cfg.bin_dir);          // 数据源目录
    gets({"MANIFEST", "MANIFEST_TSV"}, cfg.manifest);
    gets({"OUT_ROOT", "OUT", "RESULT_ROOT"}, cfg.out_root);

    geti({"BINS", "N_BINS"}, cfg.bins);
    getu64({"SEED"}, cfg.seed);

    geti({"TIER_A_N", "TIER_A"}, cfg.tier_a_n);
    geti({"TIER_B_N", "TIER_B"}, cfg.tier_b_n);

    geti({"A_ITERS"}, cfg.a_iters);
    geti({"A_BURN"}, cfg.a_burn);
    geti({"A_THIN"}, cfg.a_thin);
    geti({"A_CHAINS"}, cfg.a_chains);
    geti({"B_ITERS"}, cfg.b_iters);
    geti({"B_BURN"}, cfg.b_burn);
    geti({"B_THIN"}, cfg.b_thin);
    geti({"B_CHAINS"}, cfg.b_chains);
    geti({"C_ITERS"}, cfg.c_iters);
    geti({"C_BURN"}, cfg.c_burn);
    geti({"C_THIN"}, cfg.c_thin);
    geti({"C_CHAINS"}, cfg.c_chains);

    getf({"TAU0"}, cfg.tau0);
    getf({"K10", "KAPPA10"}, cfg.k10);
    getf({"K20", "KAPPA20"}, cfg.k20);
    getf({"LAM0", "LAMBDA0"}, cfg.lam0);
    geti({"WISHART_NU0", "NU0"}, cfg.wish_nu0);
    getf({"WISHART_B0", "B0"}, cfg.wish_b0);

    geti({"K_AUX", "M_AUX"}, cfg.k_aux);
    geti({"K_MAX_COMP", "K_MAX"}, cfg.k_max_comp);
    getu64({"M_SAMPLES"}, cfg.m_samples);

    getb({"CPU_ONLY", "CPU"}, cfg.cpu_only);
    geti({"MAX_CTX"}, cfg.max_ctx);
    getb({"SKIP_FIT"}, cfg.skip_fit);
    getb({"SKIP_SAMPLE"}, cfg.skip_sample);

    if (cfg.bins != kBins) {
        err = fmt("BINS must be %d (got %d)", kBins, cfg.bins);
        return false;
    }
    if (cfg.k_max_comp < 1 || cfg.k_max_comp > kKMax) {
        err = fmt("K_MAX_COMP must be within [1,%d]", kKMax);
        return false;
    }

    // 相对路径一律按 config 文件所在目录解析(而不是进程 CWD):
    // 这样整棵 cpp_round 目录可以整体搬盘/改名/打包,不必改任何一行;
    // 从任何目录启动 exe 也不会走丢。
    {
        std::error_code ec;
        fs::path cf = fs::absolute(fs::path(path), ec);
        fs::path base = cf.parent_path();
        auto fix = [&](std::string& s) {
            if (s.empty()) return;
            fs::path p(s);
            if (p.is_relative()) s = (base / p).lexically_normal().string();
        };
        fix(cfg.bin_dir);
        fix(cfg.manifest);
        fix(cfg.out_root);
    }
    return true;
}

}  // namespace vd
