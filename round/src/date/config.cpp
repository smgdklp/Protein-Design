#include "config.h"

#include "util.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <stdexcept>

namespace {

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
                    // 关键：不是合法 JSON 转义时保留反斜杠本身
                    out.push_back('\\');
                    out.push_back(e);
                    break;
            }
        }
        return out;
    }

    JsonValue object() {
        ++i_;  // '{'
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
        ++i_;  // '['
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

const JsonValue* findAny(const JsonValue& root, std::initializer_list<const char*> keys) {
    for (const char* k : keys) {
        const JsonValue* v = root.find(k);
        if (v != nullptr) return v;
    }
    return nullptr;
}

std::string requireString(const JsonValue& root, std::initializer_list<const char*> keys) {
    const JsonValue* v = findAny(root, keys);
    if (v == nullptr) {
        throw std::runtime_error(std::string("missing config key: ") + *keys.begin());
    }
    if (v->type != JsonValue::Type::String || v->string.empty()) {
        throw std::runtime_error(std::string("config key is not a non-empty string: ") + *keys.begin());
    }
    return v->string;
}

int asInt(const JsonValue& v, const char* key) {
    if (v.type == JsonValue::Type::Number) {
        return static_cast<int>(std::llround(v.number));
    }
    if (v.type == JsonValue::Type::String) {
        char* end = nullptr;
        long n = std::strtol(v.string.c_str(), &end, 10);
        if (end == v.string.c_str() || *end != '\0') {
            throw std::runtime_error(std::string("config key is not an integer: ") + key);
        }
        return static_cast<int>(n);
    }
    throw std::runtime_error(std::string("config key is not a number: ") + key);
}

std::filesystem::path resolvePath(const std::filesystem::path& base, const std::string& value) {
    std::filesystem::path p = utf8ToPath(value);
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

std::string jsonEscape(const std::string& s) {
    std::string out;
    out.reserve(s.size() + 8);
    for (unsigned char c : s) {
        switch (c) {
            case '"':
                out += "\\\"";
                break;
            case '\\':
                out += "\\\\";
                break;
            case '\n':
                out += "\\n";
                break;
            case '\r':
                out += "\\r";
                break;
            case '\t':
                out += "\\t";
                break;
            default:
                if (c < 0x20) {
                    char buf[8];
                    std::snprintf(buf, sizeof(buf), "\\u%04x", static_cast<unsigned>(c));
                    out += buf;
                } else {
                    out.push_back(static_cast<char>(c));
                }
                break;
        }
    }
    return out;
}

std::string jsonNumber(double v, int decimals) {
    if (std::fabs(v) < 5e-7) v = 0.0;
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.*f", decimals, v);
    return buf;
}

Config loadConfig(const std::filesystem::path& file) {
    std::ifstream in(file, std::ios::binary);
    if (!in) throw std::runtime_error("cannot open config file: " + pathToUtf8(file));
    std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    if (text.size() >= 3 && static_cast<unsigned char>(text[0]) == 0xEF &&
        static_cast<unsigned char>(text[1]) == 0xBB && static_cast<unsigned char>(text[2]) == 0xBF) {
        text.erase(0, 3);
    }

    JsonValue root = parseJson(text, "pbd_config.json");
    if (root.type != JsonValue::Type::Object) {
        throw std::runtime_error("pbd_config.json: top level must be an object");
    }

    Config cfg;
    cfg.sourceFile = file;
    const std::filesystem::path base = file.parent_path();

    cfg.pbdPath = resolvePath(base, requireString(root, {"PBD_PATH", "PDB_PATH"}));
    cfg.binPath = resolvePath(base, requireString(root, {"BIN_PATH"}));

    if (const JsonValue* v = root.find("NUM")) {
        cfg.threads = asInt(*v, "NUM");
    }
    if (cfg.threads < 1) cfg.threads = 1;
    if (cfg.threads > 64) cfg.threads = 64;

    if (const JsonValue* v = root.find("BREAK_DIST")) {
        if (v->type == JsonValue::Type::Number) {
            cfg.breakDist = v->number;
        } else if (v->type == JsonValue::Type::String) {
            cfg.breakDist = std::atof(v->string.c_str());
        } else {
            throw std::runtime_error("BREAK_DIST must be a number");
        }
    }
    if (cfg.breakDist < 0.0) cfg.breakDist = 0.0;

    if (const JsonValue* v = root.find("MIN_LEN")) {
        cfg.minLen = asInt(*v, "MIN_LEN");
        if (cfg.minLen < 1) cfg.minLen = 1;
    }

    if (const JsonValue* v = root.find("RECURSIVE")) {
        if (v->type == JsonValue::Type::Bool) {
            cfg.recursive = v->boolean;
        } else if (v->type == JsonValue::Type::Number) {
            cfg.recursive = v->number != 0.0;
        }
    }

    return cfg;
}
