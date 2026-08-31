#include "classconfig.h"

namespace SA {
    inline constexpr double Cb1 = 0.1355;
    inline constexpr double Cb2 = 0.622;
    inline constexpr double Cw1 = 3.2391;
    inline constexpr double Cw2 = 0.3;
    inline constexpr double Cw3 = 2.0;
    inline constexpr double C5 = 3.5;       // 可压缩修正,可能会被弃用
    inline constexpr double Cv1 = 7.1;
    inline constexpr double Ct3 = 1.2;
    inline constexpr double Ct4 = 0.5;
    inline constexpr double fv3 = 1.0;
    inline constexpr double kappa = 0.41;
    inline constexpr double inv_sigma = 1.5;
    inline constexpr double rmax = 10.0;
    inline constexpr double Pr = 0.71;
    inline constexpr double Prt = 0.9;

    void SA_diffusion(cc::cell_class& cell);
    void SA_source(cc::cell_class& cell);
}