#include "npy.h"

#include <cstdio>
#include <stdexcept>
#include <string>

namespace {

std::string makeHeader(const char* descr, size_t n0, size_t n1, size_t n2) {
    std::string dict = "{'descr': '";
    dict += descr;
    dict += "', 'fortran_order': False, 'shape': (";
    dict += std::to_string(n0);
    dict += ", ";
    dict += std::to_string(n1);
    dict += ", ";
    dict += std::to_string(n2);
    dict += "), }";

    // magic(6) + version(2) + header_len(2) + header 需要补齐到 64 字节整数倍
    const size_t base = 10;
    size_t total = base + dict.size() + 1;  // +1 为结尾换行
    const size_t pad = (64 - (total % 64)) % 64;
    dict.append(pad, ' ');
    dict.push_back('\n');
    return dict;
}

template <typename T>
bool writeNpy(const std::filesystem::path& file, const T* data, size_t n0, size_t n1, size_t n2,
              const char* descr) {
    std::FILE* f = nullptr;
    if (_wfopen_s(&f, file.c_str(), L"wb") != 0 || f == nullptr) return false;

    const std::string header = makeHeader(descr, n0, n1, n2);
    const unsigned char magic[6] = {0x93, 'N', 'U', 'M', 'P', 'Y'};
    const unsigned char version[2] = {1, 0};
    const uint16_t headerLen = static_cast<uint16_t>(header.size());

    bool ok = std::fwrite(magic, 1, sizeof(magic), f) == sizeof(magic) &&
              std::fwrite(version, 1, sizeof(version), f) == sizeof(version) &&
              std::fwrite(&headerLen, sizeof(headerLen), 1, f) == 1 &&
              std::fwrite(header.data(), 1, header.size(), f) == header.size();

    const size_t count = n0 * n1 * n2;
    if (ok && count != 0) {
        ok = std::fwrite(data, sizeof(T), count, f) == count;
    }
    if (std::fclose(f) != 0) ok = false;
    return ok;
}

}  // namespace

bool writeNpyF32(const std::filesystem::path& file, const float* data, size_t n0, size_t n1, size_t n2) {
    return writeNpy<float>(file, data, n0, n1, n2, "<f4");
}

bool writeNpyU32(const std::filesystem::path& file, const uint32_t* data, size_t n0, size_t n1, size_t n2) {
    return writeNpy<uint32_t>(file, data, n0, n1, n2, "<u4");
}
