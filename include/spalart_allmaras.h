#pragma once
#include "mesh.h"

namespace sa {
inline constexpr double Cb1 = 0.1355;
inline constexpr double Cb2 = 0.622;
inline constexpr double Cw1 = 3.2391;
inline constexpr double Cw2 = 0.3;
inline constexpr double Cw3 = 2.0;
inline constexpr double Cv1 = 7.1;
inline constexpr double Ct3 = 1.2;
inline constexpr double Ct4 = 0.5;
inline constexpr double kappa = 0.41;
inline constexpr double inv_sigma = 1.5;
inline constexpr double rmax = 10.0;
inline constexpr double relaxation_factor = 0.5; // 湍流方程欠松弛因子
inline constexpr double Prt = 0.9;

// 面梯度就绪后，计算黏性通量与 SA 对流、扩散系数。
void prepare_face_flux(cfd::Face& face);
// 计算当前 RK 阶段的 SA 更新值，写入 nu_tilde_next。
void advance_turbulence(cfd::Cell& cell, double coefficient);
} // namespace sa
