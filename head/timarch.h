#pragma once
#include "classconfig.h"
#include <vector>

namespace fatime {
    inline double CFL = 1.0;
    inline bool USE_GLOBAL_DT = false;
}

namespace RK {
    inline std::vector<double> RK = {0.25,1.0/6,0.375,0.5,1};
}

void local_timestep(cc::cell_class& cell);