#include "segment.h"

#include "config.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <utility>

namespace {

Vec3 sub(const Vec3& a, const Vec3& b) {
    return Vec3{a.x - b.x, a.y - b.y, a.z - b.z};
}

Vec3 cross(const Vec3& a, const Vec3& b) {
    return Vec3{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

double length(const Vec3& a) {
    return std::sqrt(dot(a, a));
}

double distance(const Vec3& a, const Vec3& b) {
    return length(sub(a, b));
}

// 标准二面角，直接返回 (cos, sin)。
// 符号约定已用真实结构验证：α 螺旋 φ≈-60°、ψ≈-45°，β 折叠 φ≈-120°、ψ≈+135°。
std::pair<double, double> dihedralCosSin(const Vec3& p1, const Vec3& p2, const Vec3& p3, const Vec3& p4) {
    const Vec3 ab = sub(p2, p1);
    const Vec3 bc = sub(p3, p2);
    const Vec3 cd = sub(p4, p3);
    const Vec3 n1 = cross(ab, bc);
    const Vec3 n2 = cross(bc, cd);
    const double bcLen = length(bc);
    if (bcLen < 1e-9 || length(n1) < 1e-9 || length(n2) < 1e-9) {
        return {0.0, 0.0};
    }
    const Vec3 bcu{bc.x / bcLen, bc.y / bcLen, bc.z / bcLen};
    const Vec3 m = cross(n1, bcu);
    const double x = dot(n1, n2);
    const double y = dot(m, n2);
    const double r = std::sqrt(x * x + y * y);
    if (r < 1e-12) return {0.0, 0.0};
    // 取负号以对齐蛋白质学界惯例（与 Biopython / PyMOL 一致）
    return {x / r, -y / r};
}

char residueToChar(const std::string& name) {
    static const struct {
        const char* name;
        char code;
    } kMap[] = {
        {"ALA", 'A'}, {"ARG", 'R'}, {"ASN", 'N'}, {"ASP", 'D'}, {"CYS", 'C'},
        {"GLN", 'Q'}, {"GLU", 'E'}, {"GLY", 'G'}, {"HIS", 'H'}, {"ILE", 'I'},
        {"LEU", 'L'}, {"LYS", 'K'}, {"MET", 'M'}, {"PHE", 'F'}, {"PRO", 'P'},
        {"SER", 'S'}, {"THR", 'T'}, {"TRP", 'W'}, {"TYR", 'Y'}, {"VAL", 'V'},
    };
    for (const auto& m : kMap) {
        if (name == m.name) return m.code;
    }
    return 'X';
}

void appendFloatArray(std::string& out, const char* key, const std::vector<double>& values) {
    out += ",\"";
    out += key;
    out += "\":[";
    for (size_t i = 0; i < values.size(); ++i) {
        if (i != 0) out.push_back(',');
        out += jsonNumber(values[i], 6);
    }
    out.push_back(']');
}

}  // namespace

std::vector<Segment> buildSegments(const std::vector<Residue>& residues,
                                   double breakDist,
                                   int minLen,
                                   SegmentStats& stats) {
    stats = SegmentStats{};
    std::vector<Segment> out;
    std::vector<const Residue*> cur;

    auto flush = [&]() {
        if (cur.empty()) return;
        if (static_cast<int>(cur.size()) >= minLen) {
            Segment seg;
            seg.chain = cur.front()->chain;
            seg.startSeq = cur.front()->seq;
            seg.endSeq = cur.back()->seq;
            seg.startIcode = cur.front()->icode;
            seg.seq.reserve(cur.size());
            seg.n.reserve(cur.size());
            seg.ca.reserve(cur.size());
            seg.c.reserve(cur.size());
            for (const Residue* r : cur) {
                seg.seq.push_back(residueToChar(r->name));
                seg.n.push_back(r->n);
                seg.ca.push_back(r->ca);
                seg.c.push_back(r->c);
            }
            out.push_back(std::move(seg));
            ++stats.keptSegments;
            stats.keptResidues += cur.size();
        } else {
            ++stats.droppedSegments;
            stats.droppedResidues += cur.size();
        }
        cur.clear();
    };

    for (const Residue& r : residues) {
        const bool valid = r.hasN && r.hasCA && r.hasC && isStandardResidueName(r.name);
        if (!valid) {
            flush();
            continue;
        }

        bool split = false;
        if (!cur.empty()) {
            const Residue* prev = cur.back();
            if (prev->chain != r.chain) {
                split = true;
            } else if (r.terBefore) {
                split = true;
            } else if (breakDist > 0.0 && distance(prev->c, r.n) > breakDist) {
                split = true;
                ++stats.distanceBreaks;
            }
        }
        if (split) flush();
        cur.push_back(&r);
    }
    flush();

    return out;
}

std::string segmentToJson(const Segment& seg,
                          const std::string& segId,
                          const std::string& pdbId,
                          const std::string& srcFile) {
    const size_t n = seg.n.size();
    std::vector<double> cosPhi(n, 0.0), sinPhi(n, 0.0), cosPsi(n, 0.0), sinPsi(n, 0.0);

    for (size_t i = 1; i < n; ++i) {
        const auto cs = dihedralCosSin(seg.c[i - 1], seg.n[i], seg.ca[i], seg.c[i]);
        cosPhi[i] = cs.first;
        sinPhi[i] = cs.second;
    }
    for (size_t i = 0; i + 1 < n; ++i) {
        const auto cs = dihedralCosSin(seg.n[i], seg.ca[i], seg.c[i], seg.n[i + 1]);
        cosPsi[i] = cs.first;
        sinPsi[i] = cs.second;
    }

    std::string out;
    out.reserve(n * 48 + 320);
    out += "{\"v\":1,\"seg\":\"";
    out += jsonEscape(segId);
    out += "\",\"pdb\":\"";
    out += jsonEscape(pdbId);
    out += "\",\"chain\":\"";
    if (seg.chain != ' ') out.push_back(seg.chain);
    out += "\",\"start\":";
    out += std::to_string(seg.startSeq);
    out += ",\"end\":";
    out += std::to_string(seg.endSeq);
    out += ",\"len\":";
    out += std::to_string(static_cast<long long>(n));
    out += ",\"seq\":\"";
    out += jsonEscape(seg.seq);
    out.push_back('"');

    appendFloatArray(out, "cos_phi", cosPhi);
    appendFloatArray(out, "sin_phi", sinPhi);
    appendFloatArray(out, "cos_psi", cosPsi);
    appendFloatArray(out, "sin_psi", sinPsi);

    out += ",\"src\":\"";
    out += jsonEscape(srcFile);
    out += "\"}\n";
    return out;
}
