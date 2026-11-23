#include "config.hpp"

#include <charconv>
#include <fstream>
#include <string_view>

#include "json_lite.hpp"

namespace pbd {

namespace {

std::string_view Trim(std::string_view s) noexcept {
    while (!s.empty() && (s.front() == ' ' || s.front() == '\t' || s.front() == '\r' || s.front() == '\n')) {
        s.remove_prefix(1);
    }
    while (!s.empty() && (s.back() == ' ' || s.back() == '\t' || s.back() == '\r' || s.back() == '\n')) {
        s.remove_suffix(1);
    }
    return s;
}

} // namespace

bool Config::Load(const std::filesystem::path& file, Config& out, std::string& err) {
    out.values.clear();
    err.clear();

    std::ifstream in(file, std::ios::binary);
    if (!in) {
        err = "打不开配置文件";
        return false;
    }
    std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    if (!in.good() && !in.eof()) {
        err = "读取配置文件失败";
        return false;
    }

    json::FlatObject flat;
    if (!json::ParseFlat(text.data(), text.size(), flat, err)) {
        return false;
    }

    for (const auto& kv : flat.items) {
        if (kv.first.empty()) {
            err = "配置里存在空键";
            return false;
        }
        const std::string_view value = Trim(kv.second);
        out.values.insert_or_assign(hash::Fnv1a(kv.first), std::string(value));
    }
    return true;
}

const std::string* Config::Get(std::uint64_t key) const noexcept {
    const auto it = values.find(key);
    return it == values.end() ? nullptr : &it->second;
}

bool Config::GetString(std::uint64_t key, std::string& v) const {
    const std::string* p = Get(key);
    if (p == nullptr) return false;
    v = *p;
    return true;
}

bool Config::GetU32(std::uint64_t key, std::uint32_t& v) const {
    const std::string* p = Get(key);
    if (p == nullptr) return false;
    const std::string_view s = Trim(*p);
    if (s.empty()) return false;
    std::uint32_t tmp = 0;
    const auto r = std::from_chars(s.data(), s.data() + s.size(), tmp);
    if (r.ec != std::errc{} || r.ptr != s.data() + s.size()) return false;
    v = tmp;
    return true;
}

} // namespace pbd
