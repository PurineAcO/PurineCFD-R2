#pragma once
#include "flow.h"

// 由温度 T（K）计算声速 a = √(γRT)，单位 m/s。
double sound_speed(double T);
namespace sutherland {

inline constexpr double Ts = 110.4;
inline constexpr double T0 = 273.15;
inline constexpr double mu0 = 1.716e-5;
// 由温度 T（K）计算 Sutherland 动力黏度，单位 Pa·s。
double dynamic_viscosity(double T);
} // namespace sutherland
