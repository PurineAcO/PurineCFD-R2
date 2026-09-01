#include "classconfig.h"

namespace jst {
    inline constexpr double k2 = 0.5;
    inline constexpr double k4 = 0.01;

    // JST用的非结构的耗散项·初始化
    void JST_dissipation_INIT(cc::cell_class& cell);

    // 更新激波探测器
    void shockwave_recognize(cc::cell_class& cell);

    // 伪Laplace算子
    void laplace_dissipation(cc::cell_class& cell);

    // JST用的非结构的耗散项
    void JST_dissipation(cc::cell_class &cell);

}

