#pragma once

#include <cstdint>
#include <filesystem>
#include <string>

// .npy v1.0（小端、C 顺序、3 维）。
bool writeNpyF32(const std::filesystem::path& file, const float* data, size_t n0, size_t n1, size_t n2);

bool writeNpyU32(const std::filesystem::path& file, const uint32_t* data, size_t n0, size_t n1, size_t n2);
