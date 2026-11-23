#pragma once
// UTF-8 字符串 <-> std::filesystem::path 的唯一转换口。
// Windows 上 path 的窄字符构造按本地代码页解释, 中文路径会烂;
// 这里统一走 u8string, 语义固定为 UTF-8。
#include <filesystem>
#include <string>
#include <string_view>

namespace pbd {

inline std::filesystem::path PathFromUtf8(std::string_view s) {
    return std::filesystem::path(
        std::u8string(reinterpret_cast<const char8_t*>(s.data()), s.size()));
}

inline std::string Utf8FromPath(const std::filesystem::path& p) {
    const std::u8string u = p.u8string();
    return std::string(reinterpret_cast<const char*>(u.data()), u.size());
}

} // namespace pbd
