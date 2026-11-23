#pragma once
// 统一键哈希: FNV-1a 64。
// 全文约定: 所有配置键、配置值、枚举取值对外一律是字符串, 一进程序就折成 uint64 哈希码,
// 之后内部只比较哈希码, 不再拿字符串做分支。
//
// 契约(跨语言/跨模块必须一致, 与 tool/src/cpp_pdbtrans 同一算法):
//   输入  : 键字符串的 UTF-8 字节, 区分大小写, 不做任何归一化/去空白
//   offset: 0xCBF29CE484222325
//   prime : 0x100000001B3
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <string>
#include <string_view>

namespace pbd::hash {

inline constexpr std::uint64_t kFnvOffset = 14695981039346656037ull; // 0xCBF29CE484222325
inline constexpr std::uint64_t kFnvPrime  = 1099511628211ull;        // 0x100000001B3

constexpr std::uint64_t Fnv1a(std::string_view s, std::uint64_t seed = kFnvOffset) noexcept {
    std::uint64_t h = seed;
    for (char ch : s) {
        h ^= static_cast<std::uint64_t>(static_cast<unsigned char>(ch));
        h *= kFnvPrime;
    }
    return h;
}

// 嵌套键: 用父键哈希做种子继续滚动, 不需要真的拼字符串
constexpr std::uint64_t Child(std::uint64_t parent, std::string_view name) noexcept {
    return Fnv1a(name, parent * kFnvPrime);
}

inline constexpr std::uint64_t operator""_h(const char* s, std::size_t n) noexcept {
    return Fnv1a(std::string_view{s, n});
}

inline std::string Hex(std::uint64_t v) {
    char buf[32];
    std::snprintf(buf, sizeof buf, "0x%016llX", static_cast<unsigned long long>(v));
    return std::string(buf);
}

} // namespace pbd::hash
