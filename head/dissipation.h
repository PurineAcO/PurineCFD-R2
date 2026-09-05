#pragma once
#include "classconfig.h"

/*说明:在一些书籍中,要求对JST的dissipation加入几何权theta_IJ,本求解器未做处理0
由于JST不会是我们的最终选择,因此这只会是一个过渡.
*/

namespace jst {
    inline constexpr double k2 = 1.0;       // 二阶阻尼
    inline constexpr double k4 = 0.05;      // 四阶阻尼

    // 激波检测器
    void shockwave_recognize(cc::cell_class& cell);
    // 本家网格的Laplace算子
    void laplace_dissipation(cc::cell_class& cell);
    // 形成JST耗散项
    void JST_dissipation(cc::cell_class &cell);
}
