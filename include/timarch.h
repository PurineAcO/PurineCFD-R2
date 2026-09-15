#pragma once
#include "classconfig.h"

// 伪时间定常
namespace fatime {
    inline double CFL = 1.0;
}

// Runge-Kutta显式时间推进
namespace RK {
    inline constexpr double RK[5] = {0.25,1.0/6,0.375,0.5,1.0};
}

// 当地时间步长, 取对流、扩散与SA源项的稳定性界
void local_timestep(cc::cell_class& cell);
