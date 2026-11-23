#pragma once

#include <filesystem>
#include <string>
#include <vector>

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct Residue {
    char chain = ' ';
    int seq = 0;
    char icode = ' ';
    std::string name;   // 三字母名；MSE 归一为 MET
    bool isMse = false;
    bool terBefore = false;  // 本残基之前出现过 TER
    bool hasN = false;
    bool hasCA = false;
    bool hasC = false;
    Vec3 n{};
    Vec3 ca{};
    Vec3 c{};
};

struct PdbReadStats {
    size_t residues = 0;
    size_t mseResidues = 0;
    size_t nonStandardResidues = 0;
};

bool isStandardResidueName(const std::string& name);

// 只读第一个 MODEL（无 MODEL 记录则读全文）。
// 原子只关心主链 N / CA / C，altLoc 取优先级最高者（空 > A > 其它）。
bool readPdbFirstModel(const std::filesystem::path& file,
                       std::vector<Residue>& out,
                       PdbReadStats& stats);
