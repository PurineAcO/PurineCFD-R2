#pragma once

#include "config.hpp"
#include <cmath>

// 形成声速
double get_sonic_velocity(double T);
// 形成能量
double get_energy(cc::physics phy);

namespace sutherland {

    inline constexpr double Ts  = 110.4;
    inline constexpr double T0  = 273.15;
    inline constexpr double mu0 = 1.716e-5;
    // Sutherland粘度
    double dynamic_viscosity(double T);
}

inline double get_sonic_velocity(double T){
    return sqrt(cc::gamma*cc::R*T);
}

inline double get_energy(cc::physics phy){
    return cc::Cv*phy.T + 0.5*(phy.u*phy.u + phy.v*phy.v);
}

namespace sutherland {
    double dynamic_viscosity(double T){
        return mu0*(T/T0)*std::sqrt(T/T0)*(T0 + Ts)/(T + Ts);
    }
}
