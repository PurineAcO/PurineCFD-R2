#pragma once
#include "mesh.h"

namespace jst {
inline constexpr double k2 = 0.5;      // 二阶阻尼
inline constexpr double k4 = 1.0 / 64; // 四阶阻尼

// 激波检测器
void compute_pressure_sensor(cfd::Cell& cell);
// 邻居与本单元守恒量差的和：Lᵢ = Σⱼ(Qⱼ - Qᵢ)
void compute_state_laplacian(cfd::Cell& cell);
// 形成JST耗散项
void compute_dissipation(cfd::Cell& cell);
} // namespace jst
