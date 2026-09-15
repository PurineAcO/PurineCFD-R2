#pragma once

namespace res {
    // 报告各物理量的全场最大变化量及位置
    void report_update(int step);
    // 最大归一化变化量(收敛判据)
    double relative_update();
}
