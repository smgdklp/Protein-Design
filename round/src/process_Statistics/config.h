#pragma once

#include <filesystem>
#include <string>
#include <utility>
#include <vector>

// 极简宽松 JSON：允许尾逗号、缺逗号、注释；字符串里只有 \" \\ \/ \uXXXX
// 会被转义，其它反斜杠原样保留（Windows 路径友好）。
struct JsonValue {
    enum class Type { Null, Bool, Number, String, Array, Object };

    Type type = Type::Null;
    bool boolean = false;
    double number = 0.0;
    std::string string;
    std::vector<JsonValue> array;
    std::vector<std::pair<std::string, JsonValue>> object;

    const JsonValue* find(const std::string& key) const;
};

JsonValue parseJson(const std::string& text, const char* what = "json");

struct Config {
    std::filesystem::path binPath;     // 输入：segments.jsonl 所在目录
    std::filesystem::path resultPath;  // 输出：方阵与摘要
    std::filesystem::path sourceFile;
};

Config loadConfig(const std::filesystem::path& file);
