#pragma once

#include <filesystem>
#include <string>

std::filesystem::path utf8ToPath(const std::string& utf8);
std::string pathToUtf8(const std::filesystem::path& p);
std::filesystem::path exeDir();
std::string nowIso8601();
