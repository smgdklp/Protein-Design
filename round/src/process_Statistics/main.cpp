#include "config.h"
#include "npy.h"

#include <windows.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr size_t kAlphabetSize = 20;
constexpr size_t kTriplets = kAlphabetSize * kAlphabetSize * kAlphabetSize;  // 8000
constexpr size_t kBins = 72;
constexpr size_t kCells = kBins * kBins;  // 5184
constexpr double kPi = 3.14159265358979323846;
constexpr double kBinWidthDeg = 5.0;

// 字母表顺序即矩阵索引顺序：idx = a*400 + b*20 + c
const char kAlphabet[] = "ACDEFGHIKLMNPQRSTVWY";

int gCode[256];

struct Options {
    std::string configArg;
    bool help = false;
    bool writeCounts = false;
    int inspectA = -1;
    int inspectB = -1;
    int inspectC = -1;
};

struct Summary {
    size_t lines = 0;
    size_t segments = 0;
    size_t residues = 0;
    size_t expectedSamples = 0;
    size_t countedSamples = 0;
    size_t badLines = 0;
    size_t badSequences = 0;
    size_t emptyMatrices = 0;
    size_t nonEmptyMatrices = 0;
    size_t normalizationViolations = 0;
    double maxCellProbability = 0.0;
    size_t medianTripletSamples = 0;
    size_t maxTripletSamples = 0;
};

fs::path exeDir() {
    std::vector<wchar_t> buf(32768);
    DWORD n = GetModuleFileNameW(nullptr, buf.data(), static_cast<DWORD>(buf.size()));
    if (n == 0 || static_cast<size_t>(n) >= buf.size()) throw std::runtime_error("GetModuleFileNameW failed");
    return fs::path(std::wstring(buf.data(), static_cast<size_t>(n))).parent_path();
}

void initCodes() {
    for (int i = 0; i < 256; ++i) gCode[i] = -1;
    for (size_t i = 0; i < kAlphabetSize; ++i) {
        gCode[static_cast<unsigned char>(kAlphabet[i])] = static_cast<int>(i);
    }
}

// 左边闭右边开：[ -180 + 5k, -180 + 5(k+1) )
// 恰好落在 +180 的点回绕到 0 号箱；-180 也落在 0 号箱。
int binOfDegrees(double angleDeg) {
    double t = (angleDeg + 180.0) / kBinWidthDeg;
    int b = static_cast<int>(std::floor(t));
    if (b < 0) b += static_cast<int>(kBins);
    if (b >= static_cast<int>(kBins)) b -= static_cast<int>(kBins);
    if (b < 0 || b >= static_cast<int>(kBins)) b = 0;
    return b;
}

bool readString(const std::string& s, size_t& pos, std::string& out) {
    if (pos >= s.size() || s[pos] != '"') return false;
    ++pos;
    out.clear();
    while (pos < s.size()) {
        const char c = s[pos++];
        if (c == '"') return true;
        if (c == '\\') return false;  // 序列里不会有转义
        out.push_back(c);
    }
    return false;
}

bool readDoubleArray(const std::string& s, size_t& pos, std::vector<double>& out) {
    if (pos >= s.size() || s[pos] != '[') return false;
    ++pos;
    out.clear();
    while (pos < s.size()) {
        const char c = s[pos];
        if (c == ']') {
            ++pos;
            return true;
        }
        if (c == ' ' || c == ',' || c == '\t') {
            ++pos;
            continue;
        }
        char* end = nullptr;
        const double v = std::strtod(s.c_str() + pos, &end);
        if (end == s.c_str() + pos) return false;
        out.push_back(v);
        pos = static_cast<size_t>(end - s.c_str());
    }
    return false;
}

bool parseSegmentLine(const std::string& line, std::string& seq, std::vector<double>& cosPhi,
                      std::vector<double>& sinPhi, std::vector<double>& cosPsi,
                      std::vector<double>& sinPsi, std::string& err) {
    if (line.find("\"v\":1") == std::string::npos) {
        err = "missing version";
        return false;
    }
    size_t pos = line.find("\"seq\":\"");
    if (pos == std::string::npos) {
        err = "missing seq";
        return false;
    }
    pos += 6;  // "seq":  之后紧跟值的开引号
    if (!readString(line, pos, seq)) {
        err = "bad seq";
        return false;
    }

    struct Field {
        const char* key;
        std::vector<double>* dst;
    };
    const Field fields[] = {
        {"\"cos_phi\":", &cosPhi},
        {"\"sin_phi\":", &sinPhi},
        {"\"cos_psi\":", &cosPsi},
        {"\"sin_psi\":", &sinPsi},
    };
    for (const Field& f : fields) {
        const size_t p = line.find(f.key, pos);
        if (p == std::string::npos) {
            err = std::string("missing ") + f.key;
            return false;
        }
        pos = p + std::strlen(f.key);
        if (!readDoubleArray(line, pos, *f.dst)) {
            err = std::string("bad array ") + f.key;
            return false;
        }
    }
    return true;
}

void logLine(std::ofstream& log, const std::string& msg) {
    const std::string line = msg + "\n";
    if (log) {
        log << line;
        log.flush();
    }
    std::cerr << line;
}

std::string nowLocal() {
    std::time_t t = std::time(nullptr);
    std::tm tm{};
    localtime_s(&tm, &t);
    char buf[32];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%S", &tm);
    return buf;
}

int parseTripletIndex(const std::string& s) {
    if (s.size() == 1) {
        const int c = gCode[static_cast<unsigned char>(s[0])];
        return c;
    }
    char* end = nullptr;
    long v = std::strtol(s.c_str(), &end, 10);
    if (end == s.c_str() || *end != '\0') return -1;
    if (v < 0 || v >= static_cast<long>(kAlphabetSize)) return -1;
    return static_cast<int>(v);
}

int run(int argc, char** argv) {
    initCodes();

    Options opt;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--config") {
            if (++i >= argc) throw std::runtime_error("--config needs a path");
            opt.configArg = argv[i];
        } else if (a == "--counts") {
            opt.writeCounts = true;
        } else if (a == "--inspect") {
            if (i + 3 >= argc) throw std::runtime_error("--inspect needs three letters (m-1 m m+1)");
            const int a1 = parseTripletIndex(argv[++i]);
            const int b1 = parseTripletIndex(argv[++i]);
            const int c1 = parseTripletIndex(argv[++i]);
            if (a1 < 0 || b1 < 0 || c1 < 0) throw std::runtime_error("--inspect: bad residue letter");
            opt.inspectA = a1;
            opt.inspectB = b1;
            opt.inspectC = c1;
        } else if (a == "--help" || a == "-h") {
            opt.help = true;
        } else {
            throw std::runtime_error("unknown argument: " + a);
        }
    }

    if (opt.help) {
        std::printf(
            "statgen - 统计 P(phi, psi | m-1, m, m+1) 的 72x72 离散分布\n"
            "\n"
            "用法: statgen [--config <config.json>] [--counts] [--inspect A B C]\n"
            "\n"
            "  --config   配置路径；缺省为 exe 同目录 config.json\n"
            "  --counts   额外输出 triplet_counts.npy（uint32 计数）\n"
            "  --inspect  打印某个三元组的 72x72 归一矩阵（CSV）\n"
            "\n"
            "输出: <RESULT_PATH>/triplet_probs.npy  shape=(8000,72,72) float32\n"
            "      残基顺序: ACDEFGHIKLMNPQRSTVWY, idx = a*400 + b*20 + c\n");
        return 0;
    }

    fs::path cfgPath;
    if (!opt.configArg.empty()) {
        cfgPath = fs::path(opt.configArg);
    } else {
        cfgPath = exeDir() / "config.json";
        if (!fs::exists(cfgPath)) cfgPath = fs::current_path() / "config.json";
    }

    const Config cfg = loadConfig(cfgPath);
    std::error_code ec;
    if (!fs::is_directory(cfg.binPath, ec)) {
        throw std::runtime_error("BIN_PATH is not a directory: " + cfg.binPath.string());
    }
    fs::create_directories(cfg.resultPath, ec);

    std::ofstream log(cfg.resultPath / "statgen.log", std::ios::binary | std::ios::app);
    logLine(log, nowLocal() + " config=" + cfgPath.string());
    logLine(log, nowLocal() + " input=" + cfg.binPath.string());
    logLine(log, nowLocal() + " output=" + cfg.resultPath.string());

    const fs::path segFile = cfg.binPath / "segments.jsonl";
    std::ifstream in(segFile, std::ios::binary);
    if (!in) throw std::runtime_error("cannot open " + segFile.string());

    // 可选的一致性提示：processed.jsonl 里 ok 的条数
    {
        std::ifstream proc(cfg.binPath / "processed.jsonl", std::ios::binary);
        size_t okCount = 0;
        std::string line;
        while (std::getline(proc, line)) {
            if (line.find("\"status\":\"ok\"") != std::string::npos) ++okCount;
        }
        if (okCount > 0) logLine(log, nowLocal() + " processed_ok=" + std::to_string(okCount));
    }

    const auto startTime = std::chrono::steady_clock::now();

    std::vector<uint32_t> counts(kTriplets * kCells, 0);
    Summary sum;

    std::string line;
    std::string seq;
    std::vector<double> cosPhi, sinPhi, cosPsi, sinPsi;
    cosPhi.reserve(4096);
    sinPhi.reserve(4096);
    cosPsi.reserve(4096);
    sinPsi.reserve(4096);

    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        ++sum.lines;

        std::string err;
        if (!parseSegmentLine(line, seq, cosPhi, sinPhi, cosPsi, sinPsi, err)) {
            ++sum.badLines;
            if (sum.badLines <= 5) logLine(log, nowLocal() + " BAD line " + std::to_string(sum.lines) + ": " + err);
            continue;
        }

        const size_t n = seq.size();
        if (cosPhi.size() != n || sinPhi.size() != n || cosPsi.size() != n || sinPsi.size() != n) {
            ++sum.badLines;
            continue;
        }
        ++sum.segments;
        sum.residues += n;
        if (n < 3) continue;

        sum.expectedSamples += n - 2;
        for (size_t i = 1; i + 1 < n; ++i) {
            const int a = gCode[static_cast<unsigned char>(seq[i - 1])];
            const int b = gCode[static_cast<unsigned char>(seq[i])];
            const int c = gCode[static_cast<unsigned char>(seq[i + 1])];
            if (a < 0 || b < 0 || c < 0) {
                ++sum.badSequences;
                break;
            }
            const double phi = std::atan2(sinPhi[i], cosPhi[i]) * 180.0 / kPi;
            const double psi = std::atan2(sinPsi[i], cosPsi[i]) * 180.0 / kPi;
            const int pi = binOfDegrees(phi);
            const int pj = binOfDegrees(psi);
            const size_t t = (static_cast<size_t>(a) * 20 + static_cast<size_t>(b)) * 20 + static_cast<size_t>(c);
            counts[t * kCells + static_cast<size_t>(pi) * kBins + static_cast<size_t>(pj)]++;
            ++sum.countedSamples;
        }

        if (sum.lines % 5000 == 0) {
            const double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime).count();
            char buf[256];
            std::snprintf(buf, sizeof(buf), "progress lines=%zu samples=%zu elapsed=%.1fs", sum.lines,
                          sum.countedSamples, el);
            logLine(log, nowLocal() + " " + buf);
        }
    }

    // 归一化
    std::vector<float> probs(kTriplets * kCells, 0.0f);
    std::vector<uint64_t> totals(kTriplets, 0);
    for (size_t t = 0; t < kTriplets; ++t) {
        const size_t base = t * kCells;
        uint64_t s = 0;
        for (size_t k = 0; k < kCells; ++k) s += counts[base + k];
        totals[t] = s;
        if (s == 0) {
            ++sum.emptyMatrices;
            continue;
        }
        ++sum.nonEmptyMatrices;
        const double inv = 1.0 / static_cast<double>(s);
        double acc = 0.0;
        for (size_t k = 0; k < kCells; ++k) {
            const float p = static_cast<float>(static_cast<double>(counts[base + k]) * inv);
            probs[base + k] = p;
            acc += p;
            if (p > sum.maxCellProbability) sum.maxCellProbability = p;
        }
        if (std::fabs(acc - 1.0) > 1e-5) ++sum.normalizationViolations;
    }

    std::vector<uint64_t> sortedTotals = totals;
    std::sort(sortedTotals.begin(), sortedTotals.end());
    sum.medianTripletSamples = static_cast<size_t>(sortedTotals[kTriplets / 2]);
    sum.maxTripletSamples = static_cast<size_t>(sortedTotals.back());

    const fs::path probsFile = cfg.resultPath / "triplet_probs.npy";
    if (!writeNpyF32(probsFile, probs.data(), kTriplets, kBins, kBins)) {
        throw std::runtime_error("failed to write " + probsFile.string());
    }
    if (opt.writeCounts) {
        const fs::path countsFile = cfg.resultPath / "triplet_counts.npy";
        if (!writeNpyU32(countsFile, counts.data(), kTriplets, kBins, kBins)) {
            throw std::runtime_error("failed to write " + countsFile.string());
        }
    }

    const double elapsed =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime).count();

    const bool ok = (sum.countedSamples == sum.expectedSamples) && (sum.badLines == 0) &&
                    (sum.badSequences == 0) && (sum.normalizationViolations == 0);

    char buf[512];
    std::snprintf(buf, sizeof(buf),
                  "DONE segments=%zu residues=%zu samples=%zu expected=%zu empty=%zu nonempty=%zu "
                  "median_per_triplet=%zu max_per_triplet=%zu max_cell_p=%.6f norm_violations=%zu "
                  "bad_lines=%zu elapsed=%.1fs result=%s",
                  sum.segments, sum.residues, sum.countedSamples, sum.expectedSamples, sum.emptyMatrices,
                  sum.nonEmptyMatrices, sum.medianTripletSamples, sum.maxTripletSamples,
                  sum.maxCellProbability, sum.normalizationViolations, sum.badLines, elapsed,
                  ok ? "OK" : "FAIL");
    logLine(log, nowLocal() + " " + buf);

    std::ofstream summaryFile(cfg.resultPath / "summary.txt", std::ios::binary | std::ios::trunc);
    if (summaryFile) {
        summaryFile << "input=" << cfg.binPath.string() << "\n"
                    << "segments=" << sum.segments << "\n"
                    << "residues=" << sum.residues << "\n"
                    << "samples=" << sum.countedSamples << "\n"
                    << "expected_samples=" << sum.expectedSamples << "\n"
                    << "matrices=8000\n"
                    << "empty_matrices=" << sum.emptyMatrices << "\n"
                    << "nonempty_matrices=" << sum.nonEmptyMatrices << "\n"
                    << "median_samples_per_triplet=" << sum.medianTripletSamples << "\n"
                    << "max_samples_per_triplet=" << sum.maxTripletSamples << "\n"
                    << "max_cell_probability=" << sum.maxCellProbability << "\n"
                    << "normalization_violations=" << sum.normalizationViolations << "\n"
                    << "bad_lines=" << sum.badLines << "\n"
                    << "bad_sequences=" << sum.badSequences << "\n"
                    << "elapsed_seconds=" << elapsed << "\n"
                    << "result=" << (ok ? "OK" : "FAIL") << "\n";
    }

    if (opt.inspectA >= 0) {
        const size_t t = static_cast<size_t>(opt.inspectA) * 400 + static_cast<size_t>(opt.inspectB) * 20 +
                         static_cast<size_t>(opt.inspectC);
        std::printf("triplet=%c%c%c idx=%zu\n", kAlphabet[opt.inspectA], kAlphabet[opt.inspectB],
                    kAlphabet[opt.inspectC], t);
        for (size_t i = 0; i < kBins; ++i) {
            for (size_t j = 0; j < kBins; ++j) {
                std::printf("%.6f", static_cast<double>(probs[t * kCells + i * kBins + j]));
                std::printf(j + 1 == kBins ? "\n" : ",");
            }
        }
    }

    return ok ? 0 : 1;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        return run(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << '\n';
        return 2;
    }
}
