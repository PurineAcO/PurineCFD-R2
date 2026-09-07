#pragma once
#include "config.h"

// 形成声速
double get_sonic_velocity(double T);
// 形成能量
double get_energy(cc::physics phy);
// 角度转换
double deg2rad(double deg);

namespace sutherland {

    inline constexpr double Ts  = 110.4;
    inline constexpr double T0  = 273.15;
    inline constexpr double mu0 = 1.716e-5;
    // Sutherland粘度
    double sutherland(double T);
}

