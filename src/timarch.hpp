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
#ifdef _OPENMP
#include <omp.h>
#endif

#define forrk(rktable) for (int z = 0; z < (int)std::size(rktable); z++)

// Runge-Kutta显式时间推进
namespace RK {
inline constexpr double RK5[5] = {0.25, 1.0 / 6, 0.375, 0.5, 1.0};
inline constexpr double RK3[3] = {0.2075, 0.5915, 1.0};
} // namespace RK

namespace CFLsup {
inline constexpr double Cv = 1.0;
inline const double Cn = (1.333 > cc::gamma / cc::Pr) ? 1.333 : cc::gamma / cc::Pr;
} // namespace CFLsup

// 当地时间步长
void local_timestep(cc::cell_class& cell);
// 稳态形式RK求解1次伪时间
void one_rans(double dt, bool urans);

inline void local_timestep(cc::cell_class& cell) {
    double lambda_c = 0.0;
    double lambda_v = 0.0;
    allface(cell) {
        // 必须用外法向速度的绝对值, 否则 lambda_c 可能为负, 双时间步会在 localdt=-2dt 处发散
        const double len = cell.faces[i]->nor.norm();
        const double un = cell.faces[i]->toface_jacobi(cell.faces[i]->phy.u, cell.faces[i]->phy.v) *
                          (2 * cell.fnorm[i] - 1) / len;
        lambda_c += (std::abs(un) + cell.faces[i]->otphy.a) * len;
        lambda_v += CFLsup::Cv * cell.faces[i]->tur.mueff / cell.faces[i]->phy.rho * CFLsup::Cn *
                    len * len / cell.vol;
    }
    cell.localdt = fatime::CFL * cell.vol / (lambda_c + lambda_v);
}

inline void one_rans(double dt, bool urans) {
#ifdef _OPENMP
    omp_set_num_threads(16);
#endif
// 复制当前伪时间步的守恒量
#pragma omp parallel for schedule(static)
    allcell icell(i).copyconver();
    forrk(RK::RK3) {
// 恢复基本物理量
#pragma omp parallel for schedule(static)
        allcell icell(i).prim();
        // 边界条件
        noslip_wall_boundary();
        far_field_boundary();
        // 更新虚拟网格
        update_ghost_field();
// 面上插值
#pragma omp parallel for schedule(static)
        allfac interpolate_mid(&iface(i));
#pragma omp parallel for schedule(static)
        allfac muscl(&iface(i));
// 建立梯度
#pragma omp parallel for schedule(static)
        allcell least_square_cell_based(icell(i));
#pragma omp parallel for schedule(static)
        allcell grad_onface(icell(i));
        // 形成对流项
        if (cc::scheme == 'J') {
#pragma omp parallel for schedule(static)
            allcell jst::shockwave_recognize(icell(i));
#pragma omp parallel for schedule(static)
            allcell jst::laplace_dissipation(icell(i));
#pragma omp parallel for schedule(static)
            allfac convect_JST(iface(i));
        } else {
#pragma omp parallel for schedule(static)
            allfac convect_ROE(iface(i));
        }
#pragma omp parallel for schedule(static)
        allcell assemble_flux(icell(i), cc::scheme);
// 形成扩散项
#pragma omp parallel for schedule(static)
        allfac SA::diffusion_SA(iface(i));
#pragma omp parallel for schedule(static)
        allcell SA::assemble_visflux(icell(i));
// 得到当地时间
#pragma omp parallel for schedule(static)
        allcell local_timestep(icell(i));
        // RK显式迭代
        if (!urans) {
#pragma omp parallel for schedule(static)
            allcell icell(i).conser =
                icell(i).conserformer - RK::RK3[z] / icell(i).vol * icell(i).localdt *
                                            (icell(i).convect - icell(i).visflux);
        } else {
#pragma omp parallel for schedule(static)
            allcell icell(i).conser =
                icell(i).conserformer -
                RK::RK3[z] / (1.0 / icell(i).localdt + 1.0 / (2 * dt)) *
                    (1.0 / icell(i).vol * (icell(i).convect - icell(i).visflux) +
                     1 / (2 * dt) * (icell(i).conser - icell(i).lastconser));
        }
    }
// 求解SA输运方程
#pragma omp parallel for schedule(static)
    allcell SA::SA_equation_after(icell(i), icell(i).localdt, dt, cc::urans);
}
