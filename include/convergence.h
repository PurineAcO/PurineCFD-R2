#pragma once

namespace convergence {
void report_update(int step);
// 各状态分量除以来流尺度后取最大更新量；这不是方程残差。
double normalized_max_update();
} // namespace convergence
