#pragma once
#include "classconfig.h"
#include <vector>

namespace fatime {
    inline double CFL = 1.0;            // CFL数
    inline constexpr double C = 1.0;    // 谱半径混合系数
}

namespace RK {
    inline std::vector<double> RK = {0.25,1.0/6,0.375,0.5,1};
}

// 计算当地时间步长
void local_timestep(cc::cell_class& cell);

// 把全局最小时间步长作为时间推进步长
void global_timestep(cc::cell_class& cell);