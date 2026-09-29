#pragma once

#include "HALO.hpp"
#include "SA.hpp"
#include "boundary.hpp"
#include "classconfig.hpp"
#include "config.hpp"
#include "convect.hpp"
#include "grad.hpp"
#include "interpolate.hpp"
#include <cmath>

#define forrk(rktable) for(int z=0;z<(int) std::size(rktable);z++)

// Runge-Kutta显式时间推进
namespace RK {
    inline constexpr double RK5[5] = {0.25,1.0/6,0.375,0.5,1.0};
    inline constexpr double RK3[3] = {0.2075,0.5915,1.0};
}

namespace CFLsup {
    inline constexpr double Cv = 1.0;
    inline const double Cn = (1.333 > cc::gamma/cc::Pr) ? 1.333 : cc::gamma/cc::Pr;
}

// 当地时间步长
void local_timestep(cc::cell_class& cell);
// 稳态形式RK求解1次伪时间
void one_rans(double dt,bool urans,cc::vec4 lastconver);

inline void local_timestep(cc::cell_class &cell){
    double lambda_c = 0.0;double lambda_v = 0.0;
    allface(cell){
        lambda_c += cell.faces[i]->toface_jacobi(cell.faces[i]->phy.u,cell.faces[i]->phy.v);
        lambda_c += cell.faces[i]->otphy.a * cell.faces[i]->nor.norm();
        lambda_v += CFLsup::Cv * cell.faces[i]->tur.mueff/cell.faces[i]->phy.rho * CFLsup::Cn * cell.faces[i]->nor.norm() *
                    cell.faces[i]->nor.norm() / cell.vol;
    }
    cell.localdt = fatime::CFL * cell.vol / (lambda_c + lambda_v);
}

inline void one_rans(double dt,bool urans){
    // 复制当前伪时间步的守恒量
    allcell icell(i).copyconver();
    forrk(RK::RK3){
        // 恢复基本物理量
        allcell icell(i).prim();
        // 边界条件
        noslip_wall_boundary();far_field_boundary();
        // 更新虚拟网格
        update_ghost_field();
        // 面上插值
        allfac interpolate_mid(&iface(i));
        allfac muscl(&iface(i));
        // 建立梯度
        allcell least_square_cell_based(icell(i));
        allcell grad_onface(icell(i));
        // 形成对流项
        allfac convect_ROE(iface(i)); 
        allcell assemble_flux(icell(i));
        // 形成扩散项
        allfac SA::diffusion_SA(iface(i));
        allcell SA::assemble_visflux(icell(i));
        // 得到当地时间
        allcell local_timestep(icell(i));
        // RK显式迭代
        if(!urans) allcell icell(i).conser = icell(i).conserformer - RK::RK3[z]/icell(i).vol*icell(i).localdt
                                            *(icell(i).convect-icell(i).visflux);
        else allcell icell(i).conser = icell(i).conserformer - RK::RK3[z]/(1.0/icell(i).localdt + 1.0/(2*dt))*
                    (1.0/icell(i).vol*(icell(i).convect-icell(i).visflux)+1/(2*dt)*(icell(i).conser-icell(i).lastconser));
    }
    // 求解SA输运方程
    allcell SA::SA_equation_after(icell(i),icell(i).localdt,dt,cc::urans);
}
