// von_dpmm / gpu.h — CUDA 宿主接口
#pragma once

#include "model.h"
#include <string>
#include <vector>

namespace vd {

struct GpuInfo {
    char   name[128] = {0};
    size_t vram = 0;
    int    cc_major = 0, cc_minor = 0;
    int    sm_count = 0;
};

bool gpu_info(GpuInfo& out, std::string& err);

// 上传 Bessel 表(宿主侧已建好)并初始化设备缓冲
bool gpu_init(const BesselTable& tab, std::string& err);
void gpu_shutdown();

// 直方图:直接统计 -> counts[kNctx*kCells]
bool gpu_histogram(const PointSet& ps, std::vector<uint32_t>& counts, std::string& err);

// DPMM-MCMC:分级拟合 + 式(7) 后验预测网格 + 压缩为参数
bool gpu_fit(const Config& cfg, const PointSet& ps, const CtxStat& st,
             std::vector<ParamRec>& recs, std::vector<float>& pred_grid, std::string& err);

// 从参数化模型做离散采样 -> counts[kNctx*kCells]
bool gpu_sample(const Config& cfg, const std::vector<ParamRec>& recs,
                std::vector<uint32_t>& counts, std::string& err);

// 参考实现(宿主,用于自检)
double host_logC(const BesselTable& tab, double k1, double k2, double lam);

// ---- CPU 后端(与 GPU 核同算法;CUDA 不可用时的退路 + 对拍基准)
bool cpu_histogram(const PointSet& ps, std::vector<uint32_t>& counts);
bool cpu_fit(const Config& cfg, const PointSet& ps, const CtxStat& st,
             std::vector<ParamRec>& recs, std::vector<float>& pred_grid, std::string& err);
bool cpu_sample(const Config& cfg, const std::vector<ParamRec>& recs,
                std::vector<uint32_t>& counts, std::string& err);

// 只用已保存的式(7) 网格重做参数压缩(不重跑 MCMC)
void model_grid(const ParamRec& R, float* out);
bool cpu_recompress(const Config& cfg, std::vector<ParamRec>& recs, const float* pred,
                    double* mean_tv, std::string& err);

}  // namespace vd

