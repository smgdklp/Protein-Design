#pragma once

#include <filesystem>
#include <string>
#include <utility>
#include <vector>

// 极简、偏宽松的 JSON 值。为本项目配置而写：
//   - 允许尾逗号、缺逗号；
//   - 允许 // 与 /* */ 注释；
//   - 字符串里只有 \" \\ \/ \uXXXX 会被转义，其它反斜杠原样保留
//     （这样 Windows 路径 "D:\round" 不会被误当成回车）。
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

std::string jsonEscape(const std::string& s);
std::string jsonNumber(double v, int decimals = 6);

struct Config {
    std::filesystem::path pbdPath;
    std::filesystem::path binPath;
    int threads = 5;
    double breakDist = 2.5;  // 相邻残基 C(i-1)-N(i) 距离阈值(A)，<=0 关闭
    int minLen = 3;          // 片段最小残基数
    bool recursive = false;  // 是否递归扫描子目录
    std::filesystem::path sourceFile;
};

Config loadConfig(const std::filesystem::path& file);
