#include "physic.h"
#include <cmath>
#include "config.h"

double get_sonic_velocity(double T){
    return sqrt(cc::gamma * cc::R * T);
}

double get_energy(cc::physics phy){
    return cc::Cv*phy.T + 0.5*(phy.u*phy.u + phy.v*phy.v);
}

double deg2rad(double deg){
    return deg/180*M_PI;
}

double sutherland::sutherland(double T){
    return sutherland::mu0 * std::pow(T/sutherland::T0, 1.5)
                          * (sutherland::T0 + sutherland::Ts)
                            / (T + sutherland::Ts);
}

