#pragma once
#include "mesh.h"

namespace jst {
inline constexpr double k2 = 0.5;      // 二阶耗散系数
inline constexpr double k4 = 1.0 / 64; // 四阶耗散系数

// 用相邻单元的压力差计算压力传感器。
void compute_pressure_sensor(cfd::Cell& cell);
// 邻居与本单元守恒量差的和：Lᵢ = Σⱼ(Qⱼ - Qᵢ)
void compute_state_laplacian(cfd::Cell& cell);
// 累加单元各面的二阶、四阶 JST 耗散通量。
void compute_dissipation(cfd::Cell& cell);
} // namespace jst
