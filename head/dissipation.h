#pragma once
#include "classconfig.h"

namespace jst {
    inline constexpr double k2 = 1.0;
    inline constexpr double k4 = 0.05;

    void shockwave_recognize(cc::cell_class& cell);
    void laplace_dissipation(cc::cell_class& cell);
    void JST_dissipation(cc::cell_class &cell);
}
