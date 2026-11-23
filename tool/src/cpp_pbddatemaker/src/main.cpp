// cpp_pbddatemaker 主程序
//   config.json 与 exe 同目录, 全部参数对外都是字符串, 进来立刻折成 FNV-1a 64 哈希码。
//   OUTPUT 按哈希分流: "SimgleAnle" -> NUM 个 SimgleAnle_maker, 每个 maker 独占一个线程。
//   过程在控制台按 UTF-8 打进度, 结束时所有池由 RAII 释放。
#include <windows.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cinttypes>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "config.hpp"
#include "hash.hpp"
#include "reader.hpp"
#include "utf8path.hpp"

namespace {

using Clock = std::chrono::steady_clock;

std::filesystem::path ExeDir() {
    std::wstring buf(32768, L'\0');
    const DWORD n = GetModuleFileNameW(nullptr, buf.data(), static_cast<DWORD>(buf.size()));
    if (n == 0) return std::filesystem::current_path();
    buf.resize(static_cast<std::size_t>(n));
    return std::filesystem::path(buf).parent_path();
}

// 单行显示用的文件名, 过长截断
std::string ShortName(const std::filesystem::path& p, std::size_t max_len) {
    std::string s = pbd::Utf8FromPath(p.filename());
    if (s.size() > max_len) {
        s.resize(max_len > 3 ? max_len - 3 : max_len);
        s += "...";
    }
    return s;
}

struct Progress {
    std::atomic<std::uint64_t> done{0};
    std::atomic<std::uint64_t> ok{0};
    std::atomic<std::uint64_t> fail{0};
    std::atomic<std::uint64_t> bytes{0};
    std::atomic<std::uint64_t> residues{0};
    std::mutex io;
};

void PrintConfig(const std::filesystem::path& config_file, const std::string& pbd_path,
                 const std::string& bin_path, std::uint32_t num, const std::string& output) {
    std::printf("pbddatemaker  |  config: %s\n", pbd::Utf8FromPath(config_file).c_str());
    std::printf("  PBD_PATH : \"%s\"\n", pbd_path.c_str());
    std::printf("  BIN_PATH : \"%s\"\n", bin_path.c_str());
    std::printf("  NUM      : \"%" PRIu32 "\"  -> %" PRIu32 " 个 worker\n", num, num);
    std::printf("  OUTPUT   : \"%s\"  (hash %s)\n", output.c_str(),
                pbd::hash::Hex(pbd::hash::Fnv1a(output)).c_str());
}

bool EnumerateInputs(const std::filesystem::path& root, const std::filesystem::path& skip_dir,
                     std::vector<std::filesystem::path>& out, std::string& err) {
    out.clear();
    err.clear();

    std::error_code ec;
    if (!std::filesystem::is_directory(root, ec)) {
        err = "PBD_PATH 不是可读目录: " + pbd::Utf8FromPath(root);
        return false;
    }

    std::filesystem::path skip = std::filesystem::weakly_canonical(skip_dir, ec);
    if (ec) {
        ec.clear();
        skip = skip_dir;
    }

    const auto opts = std::filesystem::directory_options::skip_permission_denied;
    for (std::filesystem::recursive_directory_iterator it(root, opts, ec), end; !ec && it != end;
         it.increment(ec)) {
        const std::filesystem::directory_entry& de = *it;
        std::error_code ec_local;
        if (de.is_directory(ec_local)) {
            std::error_code ec_canon;
            const std::filesystem::path canon = std::filesystem::weakly_canonical(de.path(), ec_canon);
            // 输出目录就在输入目录里, 必须剪掉, 否则会把自己吐出来的东西当输入
            if (!ec_canon && canon == skip) it.disable_recursion_pending();
            continue;
        }
        if (!de.is_regular_file(ec_local)) continue;

        const std::string name = pbd::Utf8FromPath(de.path().filename());
        if (pbd::EndsWithNoCase(name, ".pdb") || pbd::EndsWithNoCase(name, ".pdb.gz")) {
            out.push_back(de.path());
        }
    }
    if (ec) {
        err = "遍历 PBD_PATH 出错: " + ec.message();
        return false;
    }
    std::sort(out.begin(), out.end());
    return true;
}

} // namespace

int main() {
    SetConsoleOutputCP(CP_UTF8);
    std::setvbuf(stdout, nullptr, _IOFBF, 1u << 16);

    const std::filesystem::path exe_dir = ExeDir();
    const std::filesystem::path config_file = exe_dir / L"config.json";

    // ------------------------------------------------------------ 配置
    std::string err;
    pbd::Config cfg;
    if (!pbd::Config::Load(config_file, cfg, err)) {
        std::printf("[错误] 读取配置失败 (%s): %s\n", pbd::Utf8FromPath(config_file).c_str(), err.c_str());
        return 1;
    }

    std::string pbd_path;
    std::string bin_path;
    std::string output;
    std::uint32_t num = 0;

    if (!cfg.GetString(pbd::kKeyPbdPath, pbd_path)) {
        std::printf("[错误] config.json 缺少键 PBD_PATH (hash %s)\n", pbd::hash::Hex(pbd::kKeyPbdPath).c_str());
        return 1;
    }
    if (!cfg.GetString(pbd::kKeyBinPath, bin_path)) {
        std::printf("[错误] config.json 缺少键 BIN_PATH (hash %s)\n", pbd::hash::Hex(pbd::kKeyBinPath).c_str());
        return 1;
    }
    if (!cfg.GetString(pbd::kKeyOutput, output)) {
        std::printf("[错误] config.json 缺少键 OUTPUT (hash %s)\n", pbd::hash::Hex(pbd::kKeyOutput).c_str());
        return 1;
    }
    if (!cfg.GetU32(pbd::kKeyNum, num) || num == 0) {
        std::printf("[错误] NUM 必须是正整数(字符串形式), 当前读不到 (hash %s)\n",
                    pbd::hash::Hex(pbd::kKeyNum).c_str());
        return 1;
    }
    if (num > 64u) {
        std::printf("[警告] NUM=%" PRIu32 " 过大, 收到 64\n", num);
        num = 64u;
    }

    // OUTPUT 分流: 只认哈希码
    const std::uint64_t output_hash = pbd::hash::Fnv1a(output);
    if (output_hash != pbd::kOutputSimgleAnle) {
        std::printf("[错误] 不认识的 OUTPUT 取值 \"%s\" (hash %s), 目前只支持 \"SimgleAnle\"\n",
                    output.c_str(), pbd::hash::Hex(output_hash).c_str());
        return 1;
    }

    const std::filesystem::path pbd_root = pbd::PathFromUtf8(pbd_path);
    const std::filesystem::path bin_root = pbd::PathFromUtf8(bin_path);

    PrintConfig(config_file, pbd_path, bin_path, num, output);

    // ------------------------------------------------------------ 输入清单
    std::vector<std::filesystem::path> files;
    if (!EnumerateInputs(pbd_root, bin_root, files, err)) {
        std::printf("[错误] %s\n", err.c_str());
        return 1;
    }
    std::printf("  files    : %zu 个 .pdb/.pdb.gz\n", files.size());
    if (files.empty()) {
        std::printf("[警告] PBD_PATH 下没有可处理的 .pdb / .pdb.gz, 退出\n");
        return 0;
    }

    std::error_code ec;
    std::filesystem::create_directories(bin_root, ec);
    if (ec) {
        std::printf("[错误] 无法创建 BIN_PATH (%s): %s\n", pbd::Utf8FromPath(bin_root).c_str(),
                    ec.message().c_str());
        return 1;
    }
    std::printf("  output   -> %s\n\n", pbd::Utf8FromPath(bin_root).c_str());

    // ------------------------------------------------------------ 多线程
    const std::size_t total = files.size();
    std::atomic<std::size_t> next{0};
    Progress prog;

    const Clock::time_point t0 = Clock::now();

    std::vector<std::thread> workers;
    workers.reserve(num);
    for (std::uint32_t t = 0; t < num; ++t) {
        workers.emplace_back([&, t]() {
            try {
                pbd::SimgleAnle_maker maker;
                maker.Init(bin_path); // 每个 maker 只拿 BIN_PATH

                for (;;) {
                    const std::size_t idx = next.fetch_add(1, std::memory_order_relaxed);
                    if (idx >= total) break;

                    const std::filesystem::path& file = files[idx];
                    std::uint32_t rc = 0;
                    std::uint64_t ob = 0;
                    std::string werr;
                    const bool good = maker.Work(file, rc, ob, werr);

                    const std::uint64_t done_now = prog.done.fetch_add(1, std::memory_order_relaxed) + 1u;
                    const std::string label = ShortName(file, 30);

                    std::lock_guard<std::mutex> lk(prog.io);
                    if (good) {
                        prog.ok.fetch_add(1, std::memory_order_relaxed);
                        prog.bytes.fetch_add(ob, std::memory_order_relaxed);
                        prog.residues.fetch_add(rc, std::memory_order_relaxed);
                        std::printf("[%6zu/%6zu] %5.1f%%  %-30s %6" PRIu32 " res %9" PRIu64 " B  (t%" PRIu32 ")\n",
                                    idx + 1u, total, 100.0 * static_cast<double>(done_now) / static_cast<double>(total),
                                    label.c_str(), rc, ob, t);
                    } else {
                        prog.fail.fetch_add(1, std::memory_order_relaxed);
                        std::printf("[%6zu/%6zu] %5.1f%%  %-30s 失败: %s  (t%" PRIu32 ")\n",
                                    idx + 1u, total, 100.0 * static_cast<double>(done_now) / static_cast<double>(total),
                                    label.c_str(), werr.c_str(), t);
                    }
                    std::fflush(stdout);
                }
            } catch (const std::exception& e) {
                std::lock_guard<std::mutex> lk(prog.io);
                std::printf("[错误] 线程 t%" PRIu32 " 异常退出: %s\n", t, e.what());
                std::fflush(stdout);
            }
        });
    }
    for (std::thread& w : workers) w.join();

    const double secs = std::chrono::duration<double>(Clock::now() - t0).count();

    // ------------------------------------------------------------ 汇总
    const std::uint64_t ok = prog.ok.load();
    const std::uint64_t fail = prog.fail.load();
    const std::uint64_t bytes = prog.bytes.load();
    const std::uint64_t residues = prog.residues.load();
    const std::uint64_t done = prog.done.load();

    std::printf("\n完成: %" PRIu64 "/%zu 成功, %" PRIu64 " 失败, %" PRIu64 " 字节, %" PRIu64
                " 残基, 用时 %.3f s (%.1f files/s)\n",
                done, total, fail, bytes, residues, secs,
                secs > 0.0 ? static_cast<double>(done) / secs : 0.0);
    if (ok + fail != done) {
        std::printf("[警告] 计数不自洽: ok=%" PRIu64 " fail=%" PRIu64 " done=%" PRIu64 "\n", ok, fail, done);
    }
    if (fail != 0) {
        std::printf("[警告] 有 %" PRIu64 " 个文件失败, 详见上面的失败行\n", fail);
    }

    std::fflush(stdout);
    return 0;
}
