#pragma once

#include "pdb.h"

#include <string>
#include <vector>

struct Segment {
    char chain = ' ';
    int startSeq = 0;
    int endSeq = 0;
    char startIcode = ' ';
    std::string seq;
    std::vector<Vec3> n;
    std::vector<Vec3> ca;
    std::vector<Vec3> c;
};

struct SegmentStats {
    size_t keptSegments = 0;
    size_t keptResidues = 0;
    size_t droppedSegments = 0;   // 长度不足 minLen
    size_t droppedResidues = 0;
    size_t distanceBreaks = 0;    // 因 C-N 距离超阈值而切分
};

std::vector<Segment> buildSegments(const std::vector<Residue>& residues,
                                   double breakDist,
                                   int minLen,
                                   SegmentStats& stats);

// 一条 JSONL 记录（含结尾换行）。
// 约定：cos_phi[0] / sin_phi[0] 与 cos_psi[n-1] / sin_psi[n-1] 恒为 0，
// 其余位置有效（minLen >= 3 保证至少一个残基两个角都有定义）。
std::string segmentToJson(const Segment& seg,
                          const std::string& segId,
                          const std::string& pdbId,
                          const std::string& srcFile);
