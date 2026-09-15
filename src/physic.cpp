#include "physic.h"
#include <cmath>

double get_sonic_velocity(double T){
    return sqrt(cc::gamma*cc::R*T);
}

double get_energy(cc::physics phy){
    return cc::Cv*phy.T + 0.5*(phy.u*phy.u + phy.v*phy.v);
}

namespace sutherland {
    double dynamic_viscosity(double T){
        return mu0*(T/T0)*std::sqrt(T/T0)*(T0 + Ts)/(T + Ts);
    }
}
