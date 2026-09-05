#pragma once
#include "classconfig.h"
#include <vector>

// 伪时间定常
namespace fatime {
    inline double CFL = 1.0;
    inline bool USE_GLOBAL_DT = false;
}

// Runge-Kutta显式时间推进
namespace RK {
    inline std::vector<double> RK = {0.25,1.0/6,0.375,0.5,1};
}

// 计算当地时间
void local_timestep(cc::cell_class& cell);