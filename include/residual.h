#pragma once

namespace res {
    // 流动方程各量相对首次检查的下降比例低于该值即认为收敛(不含湍流工作变量)
    inline constexpr double drop_target = 1e-4;
    // 报告各物理量的全场最大变化量及位置
    void report_update(int step);
    // rho,u,v,e 相对首次检查下降比例中最差的一个, 越小越收敛
    double worst_drop();
}
