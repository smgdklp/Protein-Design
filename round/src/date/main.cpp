#include "config.h"
#include "pdb.h"
#include "segment.h"
#include "util.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <mutex>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr double kPi = 3.14159265358979323846;

struct Options {
    std::string configArg;
    int limit = -1;
    int threads = 0;
    bool rebuild = false;
    bool stats = false;
    bool help = false;
};

struct FileJob {
    fs::path path;
    std::string id;
    std::string filename;
};

struct ProcessResult {
    std::string segmentsJson;
    std::string processedLine;
    size_t segments = 0;
    size_t residues = 0;
};

class Logger {
public:
    explicit Logger(const fs::path& file) { file_.open(file, std::ios::binary | std::ios::app); }

    void info(const std::string& msg) {
        std::lock_guard<std::mutex> lk(m_);
        const std::string line = nowIso8601() + " " + msg;
        if (file_) {
            file_ << line << '\n';
            file_.flush();
        }
        std::cerr << line << '\n';
    }

private:
    std::mutex m_;
    std::ofstream file_;
};

// 多线程结果按输入序号保序落盘，保证同样的输入产出逐字节一致。
class OrderedWriter {
public:
    OrderedWriter(const fs::path& segFile, const fs::path& procFile, size_t total) : total_(total) {
        seg_.open(segFile, std::ios::out | std::ios::binary | std::ios::app);
        proc_.open(procFile, std::ios::out | std::ios::binary | std::ios::app);
        if (!seg_ || !proc_) {
            throw std::runtime_error("cannot open output files: " + pathToUtf8(segFile) + " / " +
                                     pathToUtf8(procFile) + " (errno=" + std::to_string(errno) + ")");
        }
        worker_ = std::thread([this] { run(); });
    }

    ~OrderedWriter() {
        if (worker_.joinable()) finish();
    }

    OrderedWriter(const OrderedWriter&) = delete;
    OrderedWriter& operator=(const OrderedWriter&) = delete;

    void submit(size_t idx, std::string segText, std::string procLine) {
        {
            std::lock_guard<std::mutex> lk(m_);
            buf_.emplace(idx, std::make_pair(std::move(segText), std::move(procLine)));
        }
        cv_.notify_one();
    }

    void finish() {
        {
            std::lock_guard<std::mutex> lk(m_);
            stop_ = true;
        }
        cv_.notify_one();
        if (worker_.joinable()) worker_.join();
    }

    size_t written() const { return written_; }

private:
    void run() {
        std::unique_lock<std::mutex> lk(m_);
        size_t flushTick = 0;
        for (;;) {
            cv_.wait(lk, [&] { return stop_ || buf_.find(next_) != buf_.end(); });
            auto it = buf_.find(next_);
            if (it == buf_.end()) {
                if (stop_) break;
                continue;
            }
            std::string segText = std::move(it->second.first);
            std::string procLine = std::move(it->second.second);
            buf_.erase(it);
            ++next_;
            lk.unlock();
            seg_ << segText;
            proc_ << procLine;
            if (++flushTick % 64 == 0) {
                seg_.flush();
                proc_.flush();
            }
            lk.lock();
        }
        seg_.flush();
        proc_.flush();
        written_ = next_;
    }

    std::mutex m_;
    std::condition_variable cv_;
    std::map<size_t, std::pair<std::string, std::string>> buf_;
    std::ofstream seg_;
    std::ofstream proc_;
    std::thread worker_;
    size_t total_ = 0;
    size_t next_ = 0;
    size_t written_ = 0;
    bool stop_ = false;
};

std::string toUpperAscii(std::string s) {
    for (char& c : s) {
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    }
    return s;
}

bool isPdbFile(const fs::path& p) {
    std::string ext = pathToUtf8(p.extension());
    for (char& c : ext) {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    return ext == ".pdb";
}

std::unordered_map<std::string, std::string> loadProcessed(const fs::path& file) {
    std::unordered_map<std::string, std::string> m;
    std::ifstream in(file, std::ios::binary);
    if (!in) return m;
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        try {
            const JsonValue j = parseJson(line, "processed.jsonl");
            const JsonValue* p = j.find("pdb");
            const JsonValue* s = j.find("status");
            if (p != nullptr && s != nullptr && p->type == JsonValue::Type::String &&
                s->type == JsonValue::Type::String) {
                m[p->string] = s->string;
            }
        } catch (...) {
            // 崩溃可能留下半行，忽略
        }
    }
    return m;
}

std::string buildProcessedLine(const FileJob& job,
                               const PdbReadStats& rs,
                               const SegmentStats& ss,
                               unsigned long long bytes,
                               const char* status,
                               const std::string& error) {
    std::string line = "{\"pdb\":\"";
    line += jsonEscape(job.id);
    line += "\",\"status\":\"";
    line += status;
    line.push_back('"');
    if (std::string(status) == "ok") {
        line += ",\"segments\":" + std::to_string(ss.keptSegments);
        line += ",\"residues\":" + std::to_string(ss.keptResidues);
        line += ",\"dropped_segments\":" + std::to_string(ss.droppedSegments);
        line += ",\"dropped_residues\":" + std::to_string(ss.droppedResidues);
        line += ",\"dist_breaks\":" + std::to_string(ss.distanceBreaks);
        line += ",\"mse\":" + std::to_string(rs.mseResidues);
        line += ",\"nonstandard\":" + std::to_string(rs.nonStandardResidues);
        line += ",\"bytes\":" + std::to_string(bytes);
    } else {
        line += ",\"error\":\"";
        line += jsonEscape(error);
        line.push_back('"');
    }
    line += ",\"file\":\"";
    line += jsonEscape(job.filename);
    line += "\",\"time\":\"";
    line += nowIso8601();
    line += "\"}\n";
    return line;
}

ProcessResult processFile(const FileJob& job, const Config& cfg) {
    std::vector<Residue> residues;
    PdbReadStats rs;
    if (!readPdbFirstModel(job.path, residues, rs)) {
        throw std::runtime_error("cannot open pdb file");
    }

    SegmentStats ss;
    std::vector<Segment> segments = buildSegments(residues, cfg.breakDist, cfg.minLen, ss);

    ProcessResult out;
    std::unordered_set<std::string> used;
    for (const Segment& seg : segments) {
        const std::string chainId = seg.chain == ' ' ? std::string("_") : std::string(1, seg.chain);
        std::string base = job.id + "_" + chainId + "_" + std::to_string(seg.startSeq);
        if (seg.startIcode != ' ') base.push_back(seg.startIcode);
        std::string id = base;
        int n = 1;
        while (used.count(id) != 0) {
            id = base + "#" + std::to_string(++n);
        }
        used.insert(id);
        out.segmentsJson += segmentToJson(seg, id, job.id, job.filename);
    }
    out.segments = ss.keptSegments;
    out.residues = ss.keptResidues;

    std::error_code ec;
    const unsigned long long bytes = static_cast<unsigned long long>(fs::file_size(job.path, ec));
    out.processedLine = buildProcessedLine(job, rs, ss, ec ? 0ULL : bytes, "ok", "");
    return out;
}

int runStats(const fs::path& file) {
    std::ifstream in(file, std::ios::binary);
    if (!in) {
        std::cerr << "cannot open " << pathToUtf8(file) << '\n';
        return 1;
    }

    long long segments = 0;
    long long residues = 0;
    long long violations = 0;
    long long placeholderZeros = 0;
    std::array<long long, 36> phiBins{};
    std::array<long long, 36> psiBins{};

    auto binOf = [](double angle) {
        int b = static_cast<int>(std::floor((angle + kPi) / (2.0 * kPi) * 36.0));
        if (b < 0) b = 0;
        if (b > 35) b = 35;
        return b;
    };

    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;

        JsonValue j;
        try {
            j = parseJson(line, "segments.jsonl");
        } catch (...) {
            ++violations;
            continue;
        }
        const JsonValue* seq = j.find("seq");
        const JsonValue* cp = j.find("cos_phi");
        const JsonValue* sp = j.find("sin_phi");
        const JsonValue* cq = j.find("cos_psi");
        const JsonValue* sq = j.find("sin_psi");
        if (seq == nullptr || cp == nullptr || sp == nullptr || cq == nullptr || sq == nullptr ||
            seq->type != JsonValue::Type::String || cp->type != JsonValue::Type::Array ||
            sp->type != JsonValue::Type::Array || cq->type != JsonValue::Type::Array ||
            sq->type != JsonValue::Type::Array) {
            ++violations;
            continue;
        }
        const size_t n = cp->array.size();
        if (sp->array.size() != n || cq->array.size() != n || sq->array.size() != n ||
            seq->string.size() != n) {
            ++violations;
            continue;
        }

        ++segments;
        residues += static_cast<long long>(n);

        for (size_t i = 0; i < n; ++i) {
            const double c1 = cp->array[i].number;
            const double s1 = sp->array[i].number;
            const double c2 = cq->array[i].number;
            const double s2 = sq->array[i].number;

            const double r1 = c1 * c1 + s1 * s1;
            const double r2 = c2 * c2 + s2 * s2;
            if (r1 < 1e-12) {
                ++placeholderZeros;
            } else if (std::fabs(r1 - 1.0) > 1e-4) {
                ++violations;
            }
            if (r2 < 1e-12) {
                ++placeholderZeros;
            } else if (std::fabs(r2 - 1.0) > 1e-4) {
                ++violations;
            }

            if (r1 > 0.5) phiBins[binOf(std::atan2(s1, c1))]++;
            if (r2 > 0.5) psiBins[binOf(std::atan2(s2, c2))]++;
        }
    }

    std::printf("segments: %lld\n", segments);
    std::printf("residues: %lld\n", residues);
    std::printf("consistency violations: %lld\n", violations);
    std::printf("placeholder zero slots: %lld (expected 2*segments=%lld)\n", placeholderZeros,
                segments * 2);
    std::printf("\n  bin(deg)        phi        psi\n");
    for (int b = 0; b < 36; ++b) {
        const int lo = -180 + b * 10;
        std::printf("  [%4d,%4d) %10lld %10lld\n", lo, lo + 10, phiBins[b], psiBins[b]);
    }
    return 0;
}

int run(int argc, char** argv) {
    Options opt;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--config") {
            if (++i >= argc) throw std::runtime_error("--config needs a path");
            opt.configArg = argv[i];
        } else if (a == "--limit") {
            if (++i >= argc) throw std::runtime_error("--limit needs a number");
            opt.limit = std::atoi(argv[i]);
        } else if (a == "--threads") {
            if (++i >= argc) throw std::runtime_error("--threads needs a number");
            opt.threads = std::atoi(argv[i]);
        } else if (a == "--rebuild") {
            opt.rebuild = true;
        } else if (a == "--stats") {
            opt.stats = true;
        } else if (a == "--help" || a == "-h") {
            opt.help = true;
        } else {
            throw std::runtime_error("unknown argument: " + a);
        }
    }

    if (opt.help) {
        std::printf(
            "pdbextract - 提取 PDB 中由 20 种标准氨基酸组成的连续片段\n"
            "\n"
            "用法: pdbextract [--config <pbd_config.json>] [--limit N] [--threads N]\n"
            "                  [--rebuild] [--stats] [--help]\n"
            "\n"
            "  --config   指定配置；缺省依次找 exe 同目录、当前目录的 pbd_config.json\n"
            "  --limit    只处理前 N 个待办文件（冒烟测试用）\n"
            "  --threads  覆盖配置里的 NUM\n"
            "  --rebuild  清空已有产出后台账后重跑\n"
            "  --stats    统计已有 segments.jsonl 的片段/残基/二面角分布\n");
        return 0;
    }

    fs::path cfgPath;
    if (!opt.configArg.empty()) {
        cfgPath = utf8ToPath(opt.configArg);
    } else {
        const fs::path exeCfg = exeDir() / "pbd_config.json";
        if (fs::exists(exeCfg)) {
            cfgPath = exeCfg;
        } else {
            cfgPath = fs::current_path() / "pbd_config.json";
        }
    }

    Config cfg = loadConfig(cfgPath);
    if (opt.threads > 0) {
        cfg.threads = opt.threads;
    }

    std::error_code ec;
    if (!fs::is_directory(cfg.pbdPath, ec)) {
        throw std::runtime_error("PBD_PATH is not a directory: " + pathToUtf8(cfg.pbdPath));
    }
    fs::create_directories(cfg.binPath, ec);

    const fs::path segFile = cfg.binPath / "segments.jsonl";
    const fs::path procFile = cfg.binPath / "processed.jsonl";

    if (opt.stats) {
        return runStats(segFile);
    }

    Logger log(cfg.binPath / "run.log");
    log.info("config=" + pathToUtf8(cfgPath));
    log.info("input=" + pathToUtf8(cfg.pbdPath));
    log.info("output=" + pathToUtf8(cfg.binPath));
    log.info("threads=" + std::to_string(cfg.threads) + " break_dist=" + std::to_string(cfg.breakDist) +
             " min_len=" + std::to_string(cfg.minLen));

    std::vector<fs::path> files;
    if (cfg.recursive) {
        for (auto it = fs::recursive_directory_iterator(cfg.pbdPath, ec);
             !ec && it != fs::recursive_directory_iterator(); it.increment(ec)) {
            std::error_code fec;
            if (it->is_regular_file(fec) && isPdbFile(it->path())) files.push_back(it->path());
        }
    } else {
        for (auto it = fs::directory_iterator(cfg.pbdPath, ec);
             !ec && it != fs::directory_iterator(); it.increment(ec)) {
            std::error_code fec;
            if (it->is_regular_file(fec) && isPdbFile(it->path())) files.push_back(it->path());
        }
    }
    if (ec) {
        throw std::runtime_error("failed to enumerate input directory: " + ec.message());
    }
    std::sort(files.begin(), files.end(), [](const fs::path& a, const fs::path& b) {
        return pathToUtf8(a.filename()) < pathToUtf8(b.filename());
    });
    log.info("found " + std::to_string(files.size()) + " pdb files");

    std::unordered_map<std::string, std::string> processed = loadProcessed(procFile);
    if (opt.rebuild) {
        std::error_code rec;
        fs::remove(segFile, rec);
        fs::remove(procFile, rec);
        processed.clear();
        log.info("rebuild: existing outputs removed");
    }

    std::vector<FileJob> work;
    std::unordered_set<std::string> scheduled;
    size_t skipped = 0;
    size_t duplicates = 0;
    for (const fs::path& p : files) {
        const std::string id = toUpperAscii(pathToUtf8(p.stem()));
        const auto it = processed.find(id);
        if (it != processed.end() && it->second == "ok") {
            ++skipped;
            continue;
        }
        if (!scheduled.insert(id).second) {
            ++duplicates;
            continue;
        }
        work.push_back(FileJob{p, id, pathToUtf8(p.filename())});
    }
    if (opt.limit > 0 && work.size() > static_cast<size_t>(opt.limit)) {
        work.resize(static_cast<size_t>(opt.limit));
    }

    log.info("skipped(already ok)=" + std::to_string(skipped) + " duplicate_ids=" +
             std::to_string(duplicates) + " todo=" + std::to_string(work.size()));

    if (work.empty()) {
        log.info("nothing to do");
        return 0;
    }

    const auto startTime = std::chrono::steady_clock::now();
    OrderedWriter writer(segFile, procFile, work.size());

    std::atomic<size_t> cursor{0};
    std::atomic<size_t> done{0};
    std::atomic<size_t> failed{0};
    std::atomic<size_t> totalSegments{0};
    std::atomic<size_t> totalResidues{0};

    auto workerFn = [&]() {
        for (;;) {
            const size_t idx = cursor.fetch_add(1);
            if (idx >= work.size()) break;
            const FileJob& job = work[idx];

            std::string segText;
            std::string procLine;
            size_t segs = 0;
            size_t res = 0;
            try {
                ProcessResult r = processFile(job, cfg);
                segs = r.segments;
                res = r.residues;
                segText = std::move(r.segmentsJson);
                procLine = std::move(r.processedLine);
            } catch (const std::exception& e) {
                failed.fetch_add(1);
                procLine = buildProcessedLine(job, PdbReadStats{}, SegmentStats{}, 0, "failed", e.what());
                log.info("FAILED " + job.filename + ": " + e.what());
            }

            totalSegments.fetch_add(segs);
            totalResidues.fetch_add(res);
            writer.submit(idx, std::move(segText), std::move(procLine));

            const size_t d = done.fetch_add(1) + 1;
            if (d % 200 == 0 || d == work.size()) {
                const double elapsed =
                    std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime).count();
                const double rate = elapsed > 0.0 ? static_cast<double>(d) / elapsed : 0.0;
                const double eta = rate > 0.0 ? static_cast<double>(work.size() - d) / rate : 0.0;
                char buf[256];
                std::snprintf(buf, sizeof(buf),
                              "progress %zu/%zu segments=%zu residues=%zu elapsed=%.1fs eta=%.1fs",
                              d, work.size(), totalSegments.load(), totalResidues.load(), elapsed, eta);
                log.info(buf);
            }
        }
    };

    std::vector<std::thread> threads;
    const int threadCount = std::max(1, cfg.threads);
    threads.reserve(static_cast<size_t>(threadCount));
    for (int i = 0; i < threadCount; ++i) {
        threads.emplace_back(workerFn);
    }
    for (std::thread& t : threads) {
        t.join();
    }
    writer.finish();

    const double elapsed =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime).count();
    std::error_code sec;
    const unsigned long long segBytes = static_cast<unsigned long long>(fs::file_size(segFile, sec));
    const unsigned long long procBytes = static_cast<unsigned long long>(fs::file_size(procFile, sec));

    char summary[512];
    std::snprintf(summary, sizeof(summary),
                  "DONE files=%zu ok=%zu failed=%zu segments=%zu residues=%zu elapsed=%.1fs "
                  "segments.jsonl=%llu bytes processed.jsonl=%llu bytes",
                  work.size(), work.size() - failed.load(), failed.load(), totalSegments.load(),
                  totalResidues.load(), elapsed, segBytes, procBytes);
    log.info(summary);
    if (writer.written() != work.size()) {
        log.info("WARNING: writer flushed " + std::to_string(writer.written()) + " of " +
                 std::to_string(work.size()));
    }
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        return run(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << '\n';
        return 1;
    }
}
