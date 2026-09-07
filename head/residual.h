#pragma once

namespace res {
    void report_update(int step);       // 报告: 每物理量全场最大绝对值残差及位置
    double current_residual_all();      // 各物理量均值平方残差
    double residual_absmax();           // 最大绝对值残差(收敛判据)
}
