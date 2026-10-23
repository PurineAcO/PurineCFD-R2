#pragma once
#include "classconfig.hpp"

// 边界条件参数: 由config.json读入, 供boundary.cpp使用

// 远场自由流(亚声速)

// inline void gust(double alpha,double freq){
//     FAR_DEFINE.u = 
// }

// 初场叠加绕原点的高斯涡
void add_seed(double amp);

inline void add_seed(double amp) {
    for (cc::cell_class& c : cc::CellList) {
        const double r = std::hypot(c.center.x, c.center.y);
        if (r < 0.5 || r > 3.0) continue;
        const double w = amp * std::exp(-r * r / 4.0);
        c.phy.u -= w * c.center.y;
        c.phy.v += w * c.center.x;
        c.form_otherphy();
        c.form_conservative();
    }
}