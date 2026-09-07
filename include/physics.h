#pragma once
#include "flow.h"

// 形成声速
double sound_speed(double T);
namespace sutherland {

inline constexpr double Ts = 110.4;
inline constexpr double T0 = 273.15;
inline constexpr double mu0 = 1.716e-5;
// Sutherland粘度
double dynamic_viscosity(double T);
} // namespace sutherland
