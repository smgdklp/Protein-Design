#include "config.h"

#include <windows.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <stdexcept>

namespace {

std::wstring utf8ToWide(const std::string& s) {
    if (s.empty()) return {};
    int n = MultiByteToWideChar(CP_UTF8, 0, s.data(), static_cast<int>(s.size()), nullptr, 0);
    if (n <= 0) throw std::runtime_error("invalid UTF-8 string");
    std::wstring w(static_cast<size_t>(n), L'\0');
    MultiByteToWideChar(CP_UTF8, 0, s.data(), static_cast<int>(s.size()), w.data(), n);
    return w;
}

class Parser {
public:
    Parser(const std::string& text, const char* what) : s_(text), what_(what) {}

    JsonValue parse() {
        JsonValue v = value();
        skip();
        if (i_ != s_.size()) fail("trailing characters");
        return v;
    }

private:
    const std::string& s_;
    const char* what_;
    size_t i_ = 0;

    [[noreturn]] void fail(const char* msg) {
        throw std::runtime_error(std::string(what_) + ": " + msg + " at offset " + std::to_string(i_));
    }

    char peek() const { return i_ < s_.size() ? s_[i_] : '\0'; }

    void skip() {
        for (;;) {
            while (i_ < s_.size()) {
                char c = s_[i_];
                if (c == ' ' || c == '\t' || c == '\r' || c == '\n') {
                    ++i_;
                } else {
                    break;
                }
            }
            if (i_ + 1 < s_.size() && s_[i_] == '/' && s_[i_ + 1] == '/') {
                i_ += 2;
                while (i_ < s_.size() && s_[i_] != '\n') ++i_;
                continue;
            }
            if (i_ + 1 < s_.size() && s_[i_] == '/' && s_[i_ + 1] == '*') {
                i_ += 2;
                while (i_ + 1 < s_.size() && !(s_[i_] == '*' && s_[i_ + 1] == '/')) ++i_;
                i_ = (i_ + 2 <= s_.size()) ? i_ + 2 : s_.size();
                continue;
            }
            break;
        }
    }

    JsonValue value() {
        skip();
        if (i_ >= s_.size()) fail("unexpected end of input");
        char c = s_[i_];
        if (c == '{') return object();
        if (c == '[') return array();
        if (c == '"') {
            JsonValue v;
            v.type = JsonValue::Type::String;
            v.string = parseString();
            return v;
        }
        if (c == 't' || c == 'f') return boolean();
        if (c == 'n') return nullValue();
        return number();
    }

    JsonValue boolean() {
        JsonValue v;
        v.type = JsonValue::Type::Bool;
        if (s_.compare(i_, 4, "true") == 0) {
            v.boolean = true;
            i_ += 4;
        } else if (s_.compare(i_, 5, "false") == 0) {
            v.boolean = false;
            i_ += 5;
        } else {
            fail("bad literal");
        }
        return v;
    }

    JsonValue nullValue() {
        if (s_.compare(i_, 4, "null") != 0) fail("bad literal");
        i_ += 4;
        return JsonValue{};
    }

    JsonValue number() {
        size_t start = i_;
        if (peek() == '-' || peek() == '+') ++i_;
        bool sawDigit = false;
        while (i_ < s_.size()) {
            char c = s_[i_];
            if ((c >= '0' && c <= '9') || c == '.' || c == 'e' || c == 'E' || c == '+' || c == '-') {
                if (c >= '0' && c <= '9') sawDigit = true;
                ++i_;
            } else {
                break;
            }
        }
        if (!sawDigit) fail("expected value");
        JsonValue v;
        v.type = JsonValue::Type::Number;
        try {
            v.number = std::stod(s_.substr(start, i_ - start));
        } catch (...) {
            v.number = 0.0;
        }
        return v;
    }

    static int hexVal(char c) {
        if (c >= '0' && c <= '9') return c - '0';
        if (c >= 'a' && c <= 'f') return c - 'a' + 10;
        if (c >= 'A' && c <= 'F') return c - 'A' + 10;
        return -1;
    }

    static void appendUtf8(std::string& out, unsigned cp) {
        if (cp < 0x80) {
            out.push_back(static_cast<char>(cp));
        } else if (cp < 0x800) {
            out.push_back(static_cast<char>(0xC0 | (cp >> 6)));
            out.push_back(static_cast<char>(0x80 | (cp & 0x3F)));
        } else if (cp < 0x10000) {
            out.push_back(static_cast<char>(0xE0 | (cp >> 12)));
            out.push_back(static_cast<char>(0x80 | ((cp >> 6) & 0x3F)));
            out.push_back(static_cast<char>(0x80 | (cp & 0x3F)));
        } else {
            out.push_back(static_cast<char>(0xF0 | (cp >> 18)));
            out.push_back(static_cast<char>(0x80 | ((cp >> 12) & 0x3F)));
            out.push_back(static_cast<char>(0x80 | ((cp >> 6) & 0x3F)));
            out.push_back(static_cast<char>(0x80 | (cp & 0x3F)));
        }
    }

    std::string parseString() {
        if (peek() != '"') fail("expected string");
        ++i_;
        std::string out;
        for (;;) {
            if (i_ >= s_.size()) fail("unterminated string");
            char c = s_[i_++];
            if (c == '"') break;
            if (c != '\\') {
                out.push_back(c);
                continue;
            }
            if (i_ >= s_.size()) fail("bad escape");
            char e = s_[i_++];
            switch (e) {
                case '"':
                    out.push_back('"');
                    break;
                case '\\':
                    out.push_back('\\');
                    break;
                case '/':
                    out.push_back('/');
                    break;
                case 'u': {
                    unsigned cp = 0;
                    bool ok = true;
                    if (i_ + 4 <= s_.size()) {
                        for (int k = 0; k < 4; ++k) {
                            int hv = hexVal(s_[i_ + k]);
                            if (hv < 0) {
                                ok = false;
                                break;
                            }
                            cp = cp * 16u + static_cast<unsigned>(hv);
                        }
                    } else {
                        ok = false;
                    }
                    if (ok) {
                        i_ += 4;
                        appendUtf8(out, cp);
                    } else {
                        out.push_back('\\');
                        out.push_back('u');
                    }
                    break;
                }
                default:
                    out.push_back('\\');
                    out.push_back(e);
                    break;
            }
        }
        return out;
    }

    JsonValue object() {
        ++i_;
        JsonValue v;
        v.type = JsonValue::Type::Object;
        for (;;) {
            skip();
            if (i_ >= s_.size()) fail("unterminated object");
            if (s_[i_] == '}') {
                ++i_;
                break;
            }
            if (s_[i_] == ',') {
                ++i_;
                continue;
            }
            if (s_[i_] != '"') fail("expected object key");
            std::string key = parseString();
            skip();
            if (peek() != ':') fail("expected ':'");
            ++i_;
            JsonValue val = value();
            v.object.emplace_back(std::move(key), std::move(val));
            skip();
            if (peek() == ',') ++i_;
        }
        return v;
    }

    JsonValue array() {
        ++i_;
        JsonValue v;
        v.type = JsonValue::Type::Array;
        for (;;) {
            skip();
            if (i_ >= s_.size()) fail("unterminated array");
            if (s_[i_] == ']') {
                ++i_;
                break;
            }
            if (s_[i_] == ',') {
                ++i_;
                continue;
            }
            v.array.push_back(value());
            skip();
            if (peek() == ',') ++i_;
        }
        return v;
    }
};

std::string requireString(const JsonValue& root, const char* key) {
    const JsonValue* v = root.find(key);
    if (v == nullptr) throw std::runtime_error(std::string("missing config key: ") + key);
    if (v->type != JsonValue::Type::String || v->string.empty()) {
        throw std::runtime_error(std::string("config key is not a non-empty string: ") + key);
    }
    return v->string;
}

std::filesystem::path resolvePath(const std::filesystem::path& base, const std::string& value) {
    std::filesystem::path p(utf8ToWide(value));
    if (p.is_relative()) p = base / p;
    return p.lexically_normal();
}

}  // namespace

const JsonValue* JsonValue::find(const std::string& key) const {
    if (type != Type::Object) return nullptr;
    for (const auto& kv : object) {
        if (kv.first == key) return &kv.second;
    }
    return nullptr;
}

JsonValue parseJson(const std::string& text, const char* what) {
    Parser p(text, what);
    return p.parse();
}

Config loadConfig(const std::filesystem::path& file) {
    std::ifstream in(file, std::ios::binary);
    if (!in) throw std::runtime_error("cannot open config file: " + file.string());
    std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    if (text.size() >= 3 && static_cast<unsigned char>(text[0]) == 0xEF &&
        static_cast<unsigned char>(text[1]) == 0xBB && static_cast<unsigned char>(text[2]) == 0xBF) {
        text.erase(0, 3);
    }

    const JsonValue root = parseJson(text, "config.json");
    if (root.type != JsonValue::Type::Object) {
        throw std::runtime_error("config.json: top level must be an object");
    }

    Config cfg;
    cfg.sourceFile = file;
    const std::filesystem::path base = file.parent_path();
    cfg.binPath = resolvePath(base, requireString(root, "BIN_PATH"));
    cfg.resultPath = resolvePath(base, requireString(root, "RESULT_PATH"));
    return cfg;
}
