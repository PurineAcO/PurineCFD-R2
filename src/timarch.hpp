#pragma once

#include "classconfig.hpp"
#include "SA.hpp"
#include "physic.hpp"
#include <cmath>

// Runge-Kutta显式时间推进
namespace RK {
    inline constexpr double RK[5] = {0.25,1.0/6,0.375,0.5,1.0};
}

// 当地时间步长, 取对流、扩散与SA源项的稳定性界
void local_timestep(cc::cell_class& cell);

inline void local_timestep(cc::cell_class& cell){
    const double dsx = cell.proj.x;
    const double dsy = cell.proj.y;
    double convective_x = (std::abs(cell.phy.u) + cell.phy.a)*dsx;
    double convective_y = (std::abs(cell.phy.v) + cell.phy.a)*dsy;
    double diffusivity_bound =
        sutherland::dynamic_viscosity(cell.phy.T)/cell.phy.rho + cell.tur.miubl;
    double viscous_bound = 4.0*diffusivity_bound*(dsx*dsx + dsy*dsy)*cell.invvol/cc::Pr;
    double sa_source_bound =
        2.0*SA::Cw1*cell.tur.miubl*(1.0/(cell.tur.sad*cell.tur.sad));
    cell.localdt = fatime::CFL/((convective_x + convective_y + viscous_bound)*cell.invvol +
                                sa_source_bound);
}
