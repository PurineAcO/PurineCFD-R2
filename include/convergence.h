#pragma once

namespace convergence {
void report_update(int step);
// 返回相邻两次检查之间的最大归一化状态更新量 max|Δφ|/φ_ref。
double normalized_max_update();
} // namespace convergence
