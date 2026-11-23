#pragma once
// 极简 JSON: 只认一层扁平对象, 键必须带引号; 值可以是字符串, 也可以是数字/true/false,
// 后者按字面文本收成字符串(配置项对外一律字符串)。嵌套对象/数组直接报错, 不静默吞。
//
// 容忍: UTF-8 BOM、注释 // ... 与 /* ... */、尾随逗号、键值之间的任意空白。
#include <cstddef>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace pbd::json {

struct FlatObject {
    std::vector<std::pair<std::string, std::string>> items; // (键, 值的字面文本)

    const std::string* Find(std::string_view key) const noexcept {
        for (const auto& kv : items) {
            if (kv.first == key) return &kv.second;
        }
        return nullptr;
    }
};

namespace detail {

inline bool IsSpace(char c) noexcept {
    return c == ' ' || c == '\t' || c == '\r' || c == '\n';
}

struct Cursor {
    const char* p;
    const char* e;

    bool AtEnd() const noexcept { return p >= e; }

    void SkipSpace() noexcept {
        for (;;) {
            while (p < e && IsSpace(*p)) ++p;
            if (p + 1 < e && p[0] == '/' && p[1] == '/') {
                p += 2;
                while (p < e && *p != '\n') ++p;
                continue;
            }
            if (p + 1 < e && p[0] == '/' && p[1] == '*') {
                p += 2;
                while (p + 1 < e && !(p[0] == '*' && p[1] == '/')) ++p;
                p = (p + 1 < e) ? p + 2 : e;
                continue;
            }
            break;
        }
    }
};

inline bool AppendUtf8(std::string& out, unsigned cp) {
    if (cp < 0x80u) {
        out.push_back(static_cast<char>(cp));
    } else if (cp < 0x800u) {
        out.push_back(static_cast<char>(0xC0u | (cp >> 6)));
        out.push_back(static_cast<char>(0x80u | (cp & 0x3Fu)));
    } else if (cp < 0x10000u) {
        out.push_back(static_cast<char>(0xE0u | (cp >> 12)));
        out.push_back(static_cast<char>(0x80u | ((cp >> 6) & 0x3Fu)));
        out.push_back(static_cast<char>(0x80u | (cp & 0x3Fu)));
    } else {
        out.push_back(static_cast<char>(0xF0u | (cp >> 18)));
        out.push_back(static_cast<char>(0x80u | ((cp >> 12) & 0x3Fu)));
        out.push_back(static_cast<char>(0x80u | ((cp >> 6) & 0x3Fu)));
        out.push_back(static_cast<char>(0x80u | (cp & 0x3Fu)));
    }
    return true;
}

inline bool ParseHex4(const char* p, const char* e, unsigned& v) noexcept {
    if (e - p < 4) return false;
    v = 0;
    for (int i = 0; i < 4; ++i) {
        const char c = p[i];
        unsigned d;
        if (c >= '0' && c <= '9')      d = static_cast<unsigned>(c - '0');
        else if (c >= 'a' && c <= 'f') d = static_cast<unsigned>(c - 'a') + 10u;
        else if (c >= 'A' && c <= 'F') d = static_cast<unsigned>(c - 'A') + 10u;
        else return false;
        v = v * 16u + d;
    }
    return true;
}

inline bool ParseString(Cursor& c, std::string& out, std::string& err) {
    if (c.AtEnd() || *c.p != '"') {
        err = "期望字符串起始引号";
        return false;
    }
    ++c.p;
    out.clear();
    while (c.p < c.e) {
        const char ch = *c.p;
        if (ch == '"') {
            ++c.p;
            return true;
        }
        if (ch != '\\') {
            if (static_cast<unsigned char>(ch) < 0x20u) {
                err = "字符串里出现未转义控制字符";
                return false;
            }
            out.push_back(ch);
            ++c.p;
            continue;
        }
        ++c.p;
        if (c.p >= c.e) {
            err = "反斜杠后被截断";
            return false;
        }
        const char esc = *c.p++;
        switch (esc) {
        case '"':  out.push_back('"');  break;
        case '\\': out.push_back('\\'); break;
        case '/':  out.push_back('/');  break;
        case 'b':  out.push_back('\b'); break;
        case 'f':  out.push_back('\f'); break;
        case 'n':  out.push_back('\n'); break;
        case 'r':  out.push_back('\r'); break;
        case 't':  out.push_back('\t'); break;
        case 'u': {
            unsigned cp = 0;
            if (!ParseHex4(c.p, c.e, cp)) {
                err = "\\u 转义不是 4 位十六进制";
                return false;
            }
            c.p += 4;
            if (cp >= 0xD800u && cp <= 0xDBFFu) { // 代理对
                if (c.e - c.p < 6 || c.p[0] != '\\' || c.p[1] != 'u') {
                    err = "高位代理后没有低位代理";
                    return false;
                }
                unsigned lo = 0;
                if (!ParseHex4(c.p + 2, c.e, lo) || lo < 0xDC00u || lo > 0xDFFFu) {
                    err = "低位代理非法";
                    return false;
                }
                c.p += 6;
                cp = 0x10000u + ((cp - 0xD800u) << 10) + (lo - 0xDC00u);
            }
            AppendUtf8(out, cp);
            break;
        }
        default:
            err = std::string("未知转义 \\") + esc;
            return false;
        }
    }
    err = "字符串没有闭合";
    return false;
}

// 数字 / true / false / null 一律按字面文本收下
inline bool ParseBare(Cursor& c, std::string& out) {
    const char* b = c.p;
    while (c.p < c.e) {
        const char ch = *c.p;
        const bool stop = ch == ',' || ch == '}' || ch == ']' || IsSpace(ch);
        if (stop) break;
        ++c.p;
    }
    if (c.p == b) return false;
    out.assign(b, static_cast<std::size_t>(c.p - b));
    return true;
}

} // namespace detail

inline bool ParseFlat(const char* data, std::size_t size, FlatObject& out, std::string& err) {
    out.items.clear();
    err.clear();

    const char* p = data;
    const char* e = data + size;
    // UTF-8 BOM
    if (size >= 3 && static_cast<unsigned char>(p[0]) == 0xEF &&
        static_cast<unsigned char>(p[1]) == 0xBB && static_cast<unsigned char>(p[2]) == 0xBF) {
        p += 3;
    }

    detail::Cursor c{p, e};
    c.SkipSpace();
    if (c.AtEnd() || *c.p != '{') {
        err = "根节点必须是 JSON 对象 {...}";
        return false;
    }
    ++c.p;

    for (;;) {
        c.SkipSpace();
        if (c.AtEnd()) {
            err = "对象没有闭合";
            return false;
        }
        if (*c.p == '}') { // 含尾随逗号
            ++c.p;
            break;
        }
        std::string key;
        if (!detail::ParseString(c, key, err)) return false;

        c.SkipSpace();
        if (c.AtEnd() || *c.p != ':') {
            err = "键 \"" + key + "\" 后缺少冒号";
            return false;
        }
        ++c.p;
        c.SkipSpace();

        std::string value;
        if (c.AtEnd()) {
            err = "键 \"" + key + "\" 缺少值";
            return false;
        }
        if (*c.p == '"') {
            if (!detail::ParseString(c, value, err)) return false;
        } else if (*c.p == '{' || *c.p == '[') {
            err = "不支持嵌套对象/数组(键 \"" + key + "\")";
            return false;
        } else if (!detail::ParseBare(c, value)) {
            err = "键 \"" + key + "\" 的值无法解析";
            return false;
        }
        out.items.emplace_back(std::move(key), std::move(value));

        c.SkipSpace();
        if (c.AtEnd()) {
            err = "对象没有闭合";
            return false;
        }
        if (*c.p == ',') {
            ++c.p;
            continue;
        }
        if (*c.p == '}') {
            ++c.p;
            break;
        }
        err = "键值后缺少逗号或右花括号";
        return false;
    }

    c.SkipSpace();
    if (!c.AtEnd()) {
        err = "根对象之后还有多余内容";
        return false;
    }
    return true;
}

} // namespace pbd::json
