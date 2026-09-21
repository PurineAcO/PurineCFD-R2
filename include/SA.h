#pragma once
#include "classconfig.h"

namespace SA {
    inline constexpr double Cb1 = 0.1355;
    inline constexpr double Cb2 = 0.622;
    inline constexpr double Cw1 = 3.2391;
    inline constexpr double Cw2 = 0.3;
    inline constexpr double Cw3 = 2.0;
    inline constexpr double Cv1 = 7.1;
    inline constexpr double Ct3 = 1.2;
    inline constexpr double Ct4 = 0.5;
    inline constexpr double kappa = 0.41;
    inline constexpr double inv_sigma = 1.5;
    inline constexpr double rmax = 10.0;
    inline constexpr double relax = 0.5;    // 湍流方程欠松弛因子
    inline constexpr double Prt = 0.9;
    inline constexpr double C5 = 3.5;

    // 形成湍流扩散项和黏性通量
    void diffusion_SA(cc::face_class& face);
    // RK子步内求解湍流方程(非守恒形式)
    void SA_equation_RK(cc::cell_class& cell, double rk);
}
