#include "pdb.h"

#include <cstdlib>
#include <fstream>
#include <string_view>
#include <unordered_map>

namespace {

const std::string_view kStandardNames[] = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
};

std::string trim(const std::string& s) {
    size_t b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return {};
    size_t e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
}

double toDouble(const std::string& s) {
    try {
        return std::stod(s);
    } catch (...) {
        return 0.0;
    }
}

int parseResidueSeq(const std::string& raw, int lastSeq) {
    char* end = nullptr;
    long v = std::strtol(raw.c_str(), &end, 10);
    if (end == raw.c_str()) return lastSeq + 1;
    return static_cast<int>(v);
}

int altPriority(char alt) {
    return alt == ' ' ? 3 : (alt == 'A' ? 2 : 1);
}

}  // namespace

bool isStandardResidueName(const std::string& name) {
    for (std::string_view n : kStandardNames) {
        if (n == name) return true;
    }
    return false;
}

bool readPdbFirstModel(const std::filesystem::path& file,
                       std::vector<Residue>& out,
                       PdbReadStats& stats) {
    std::ifstream in(file, std::ios::binary);
    if (!in) return false;

    out.clear();
    stats = PdbReadStats{};

    std::string line;
    std::string curKey;
    Residue cur;
    bool curStarted = false;
    bool pendingTer = false;
    int modelCount = 0;
    int lastSeq = 0;
    char altOf[3] = {' ', ' ', ' '};
    double occOf[3] = {-1.0, -1.0, -1.0};

    auto flushResidue = [&]() {
        if (curStarted) {
            out.push_back(cur);
            curStarted = false;
            curKey.clear();
        }
    };

    auto acceptAtom = [&](int which, char alt, double occ, const Vec3& pos) {
        bool* has = which == 0 ? &cur.hasN : (which == 1 ? &cur.hasCA : &cur.hasC);
        Vec3* dst = which == 0 ? &cur.n : (which == 1 ? &cur.ca : &cur.c);
        if (!*has) {
            *has = true;
            *dst = pos;
            altOf[which] = alt;
            occOf[which] = occ;
            return;
        }
        const int p = altPriority(alt);
        const int cp = altPriority(altOf[which]);
        if (p > cp || (p == cp && occ > occOf[which])) {
            *dst = pos;
            altOf[which] = alt;
            occOf[which] = occ;
        }
    };

    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();

        if (line.rfind("MODEL", 0) == 0) {
            ++modelCount;
            if (modelCount >= 2) break;
            continue;
        }
        if (line.rfind("ENDMDL", 0) == 0) break;

        if (line.rfind("TER", 0) == 0) {
            flushResidue();
            pendingTer = true;
            continue;
        }

        if (line.size() < 54) continue;
        const bool isAtom = line.rfind("ATOM", 0) == 0;
        const bool isHet = line.rfind("HETATM", 0) == 0;
        if (!isAtom && !isHet) continue;

        const std::string resName = trim(line.substr(17, 3));
        // 只保留标准氨基酸；HETATM 里只有 MSE（硒代甲硫氨酸）放行。
        if (isHet && resName != "MSE") continue;

        const bool mse = (resName == "MSE");
        const std::string name = mse ? std::string("MET") : resName;
        const char chain = line[21];
        const std::string seqRaw = trim(line.substr(22, 4));
        const char icode = line[26];
        const char alt = line[16];

        std::string key;
        key.reserve(16);
        key.push_back(chain);
        key.push_back('|');
        key += seqRaw;
        key.push_back('|');
        key.push_back(icode);
        key.push_back('|');
        key += name;

        if (key != curKey) {
            flushResidue();
            curKey = std::move(key);
            cur = Residue{};
            cur.chain = chain;
            cur.icode = icode;
            cur.name = name;
            cur.isMse = mse;
            cur.terBefore = pendingTer;
            pendingTer = false;
            cur.seq = parseResidueSeq(seqRaw, lastSeq);
            lastSeq = cur.seq;
            curStarted = true;
            altOf[0] = altOf[1] = altOf[2] = ' ';
            occOf[0] = occOf[1] = occOf[2] = -1.0;
        }

        const std::string atomName = trim(line.substr(12, 4));
        int which = -1;
        if (atomName == "N") {
            which = 0;
        } else if (atomName == "CA") {
            which = 1;
        } else if (atomName == "C") {
            which = 2;
        }
        if (which < 0) continue;

        Vec3 pos;
        pos.x = toDouble(line.substr(30, 8));
        pos.y = toDouble(line.substr(38, 8));
        pos.z = toDouble(line.substr(46, 8));
        const double occ = line.size() >= 60 ? toDouble(line.substr(54, 6)) : 0.0;
        acceptAtom(which, alt, occ, pos);
    }

    flushResidue();

    stats.residues = out.size();
    for (const Residue& r : out) {
        if (r.isMse) ++stats.mseResidues;
        if (!isStandardResidueName(r.name)) ++stats.nonStandardResidues;
    }
    return true;
}
