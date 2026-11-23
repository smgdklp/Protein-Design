#pragma once
// 配置读取: config.json 与 exe 同目录。
// 对外一切参数都是字符串, 进来立刻折成 uint64 哈希码(见 hash.hpp), 内部只认哈希码。
// 数值参数(NUM)额外用 from_chars 解析出真实数字, 哈希只负责身份。
#include <cstdint>
#include <filesystem>
#include <string>
#include <unordered_map>

#include "hash.hpp"

namespace pbd {

// 配置键的哈希码(编译期算好, 不存字符串)
inline constexpr std::uint64_t kKeyPbdPath = hash::Fnv1a("PBD_PATH");
inline constexpr std::uint64_t kKeyBinPath = hash::Fnv1a("BIN_PATH");
inline constexpr std::uint64_t kKeyNum     = hash::Fnv1a("NUM");
inline constexpr std::uint64_t kKeyOutput  = hash::Fnv1a("OUTPUT");

// "OUTPUT" 的合法取值
inline constexpr std::uint64_t kOutputSimgleAnle = hash::Fnv1a("SimgleAnle");

struct Config {
    std::unordered_map<std::uint64_t, std::string> values; // 键哈希 -> 值的原始文本

    // 读文件 + 解析 + 折哈希; 失败时 err 带具体原因
    static bool Load(const std::filesystem::path& file, Config& out, std::string& err);

    const std::string* Get(std::uint64_t key) const noexcept;

    bool GetString(std::uint64_t key, std::string& v) const;

    // 严格十进制无符号, 允许首尾空白(值本身已被去引号)
    bool GetU32(std::uint64_t key, std::uint32_t& v) const;
};

} // namespace pbd
