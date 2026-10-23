#pragma once

#include "config.hpp"
#include "classconfig.hpp"
#include <cmath>

// 形成声速
double get_sonic_velocity(double T);
// 形成能量
double get_energy(cc::vecp phy);
// 动压
double dynamic_pressure(cc::vecp refp);
// 壁面力系数
void wall_forces(double& cl, double& cd, cc::vecp refp);
// 检查坏单元
int bad_cells();

namespace sutherland {

inline constexpr double Ts = 110.4;
inline constexpr double T0 = 273.15;
inline constexpr double mu0 = 1.716e-5;
// Sutherland粘度
double dynamic_viscosity(double T);
} // namespace sutherland

inline double get_sonic_velocity(double T) {
    return sqrt(cc::gamma * cc::R * T);
}

inline double get_energy(cc::vecp phy) {
    return cc::Cv * phy.T + 0.5 * (phy.u * phy.u + phy.v * phy.v);
}

namespace sutherland {
inline double dynamic_viscosity(double T) {
    return mu0 * (T / T0) * std::sqrt(T / T0) * (T0 + Ts) / (T + Ts);
}
} // namespace sutherland

// 计算附加物理量e,a,p,mu
inline void cc::otphy::form_otphy(cc::vecp phy) {
    e = get_energy(phy);
    a = get_sonic_velocity(phy.T);
    p = phy.rho * cc::R * phy.T;
    mu = sutherland::dynamic_viscosity(phy.T);
}

inline double dynamic_pressure(cc::vecp refp) {
    // const double rho_inf = refp.otphy.p / (cc::R * refp.phy.T);
    return 0.5 * refp.rho * (refp.u*refp.u+refp.v*refp.v);
}

inline void wall_forces(double& cl, double& cd, cc::vecp refp) {
    const double q = dynamic_pressure(refp);
    double fx = 0.0, fy = 0.0;
    for (cc::face_class* w : cc::WallFaces) {
        const cc::cell_class* c = cc::boundary_findcell(w);
        // s*nor 由流体指向物面
        const double s = cc::dot(w->nor, w->mid - c->center) > 0 ? 1.0 : -1.0;
        fx += s * (w->otphy.p * w->nor.x - w->visflux.x);
        fy += s * (w->otphy.p * w->nor.y - w->visflux.y);
    }
    const double a = std::atan2(refp.v, refp.u);
    cl = (-fx * std::sin(a) + fy * std::cos(a)) / q;
    cd = (fx * std::cos(a) + fy * std::sin(a)) / q;
}

inline int bad_cells() {
    int nbad = 0;
    for (const cc::cell_class& c : cc::CellList) {
        const double ke =
            0.5 * (c.conser.x * c.conser.x + c.conser.y * c.conser.y) / c.conser.c;
        const double T = (c.conser.e - ke) / (c.conser.c * cc::Cv);
        const bool bad = !std::isfinite(c.conser.c) || c.conser.c <= 0.0 ||
                         !std::isfinite(c.conser.e) || !std::isfinite(c.tur.miubl) ||
                         !(T > 0.0) || !std::isfinite(T);
        if (bad) {
            if (nbad < 100)
                printf("  bad cell #%d (%.6g, %.6g)\n", c.index, c.center.x, c.center.y);
            ++nbad;
        }
    }
    if (nbad > 100) printf("  ... 共 %d 个坏单元\n", nbad);
    return nbad;
}