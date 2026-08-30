#include <cmath>
#include "config.h"

double deg2rad(double deg){return deg/180*M_PI;}

double get_sonic_velocity(double T){return sqrt(cc::gamma * cc::R * T);}

double SutherLand(double T){return cc::mu_ref * pow((T/cc::T_ref),1.5) * (cc::T_ref+cc::T_s)/(T+cc::T_s);}

double get_energy(cc::physics phy){return phy.p/(phy.rho*(cc::gamma - 1.0)) + 0.5*(phy.u*phy.u + phy.v*phy.v);}