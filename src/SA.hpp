#pragma once

#include "classconfig.hpp"
#include "config.hpp"
#include "physic.hpp"
#include <cmath>

/*
根据Fluent的求解手段,miubl将会在一个RK内迭代步以后再进行更新,而不是跟随RK
*/

namespace SA {
    inline constexpr double Cb1 = 0.1355;
    inline constexpr double Cb2 = 0.622;
    inline constexpr double Cw1 = 3.2391;
    inline constexpr double Cw2 = 0.3;
    inline constexpr double Cw3 = 2.0;
    inline constexpr double Cv1 = 7.1;
    inline constexpr double Ct3 = 1.2;
    inline constexpr double Ct4 = 0.5;
    inline constexpr double kappa = 0.41;
    inline constexpr double inv_sigma = 1.5;
    inline constexpr double rmax = 10.0;
    inline constexpr double relax = 0.8;
    inline constexpr double Prt = 0.9;
    inline constexpr double C5 = 3.5;

    // 形成粘性通量
    void diffusion_SA(cc::face_class& face);
    // 组装粘性通量
    void assemble_visflux(cc::cell_class& cell);
    // 形成源项
    double source_SA(const cc::cell_class &cell);
    // RK子步内求解湍流方程(非守恒形式)
    void SA_equation_after(cc::cell_class& cell,double dtau,double dt,bool urans);
}

static double _chi(double rho,double mu,double miubl){
    return rho*(miubl > 0.0 ? miubl : 0.0)/mu;
}

static double _fv1(double chi){
    return (chi*chi*chi)/(chi*chi*chi + SA::Cv1*SA::Cv1*SA::Cv1);
}

static double _ft2(double chi){
    return SA::Ct3*std::exp(-SA::Ct4*chi*chi);
}

static double _fv2(double chi){
    return 1 - chi/(1 + chi*_fv1(chi));
}

static double _g(double r){
    return r + SA::Cw2*(r*r*r*r*r*r - r);
}

static double _fw(double g){
    constexpr double cw3_squared = SA::Cw3*SA::Cw3;
    constexpr double cw3_sixth = cw3_squared*cw3_squared*cw3_squared;
    const double g_squared = g*g;
    const double g_sixth = g_squared*g_squared*g_squared;
    return g*std::pow((1 + cw3_sixth)/(cw3_sixth + g_sixth),1.0/6);
}

inline double SA::source_SA(const cc::cell_class& cell){
    double mu = cell.otphy.mu;
    double chi = _chi(cell.phy.rho,mu,cell.tur.miubl);
    double ft2 = _ft2(chi);
    double fv2 = _fv2(chi);
    double vorticity = std::abs(cell.phgrad.vgrad.x - cell.phgrad.ugrad.y);
    const double scaled_nu = cell.tur.miubl*(1.0/(cell.tur.sad*cell.tur.sad))/(SA::kappa*SA::kappa);
    double modified_vorticity =
        std::max(vorticity + fv2*scaled_nu,std::max(0.3*vorticity,1e-20));
    double production = SA::Cb1*(1 - ft2)*modified_vorticity*cell.phy.rho*cell.tur.miubl;
    double r = std::min(scaled_nu/modified_vorticity,SA::rmax);
    double g = _g(r);
    double destruction = cell.phy.rho*
                         (SA::Cw1*_fw(g) - SA::Cb1/SA::kappa/SA::kappa*ft2)*
                         cell.tur.miubl*cell.tur.miubl*(1.0/(cell.tur.sad*cell.tur.sad));
    double gradient_source = SA::Cb2*SA::inv_sigma*cell.phy.rho*
                             cc::dot(cell.tur.miublgrad,cell.tur.miublgrad);
    // 部分论文中引入了可压缩性修正
    double S2 = 2*cell.phgrad.ugrad.x*cell.phgrad.ugrad.x + 2*cell.phgrad.vgrad.y*cell.phgrad.vgrad.y +
                (cell.phgrad.ugrad.y + cell.phgrad.vgrad.x)*(cell.phgrad.ugrad.y + cell.phgrad.vgrad.x);
    double compressible = SA::C5 * cell.phy.rho * cell.tur.miubl * cell.tur.miubl * S2 / (cc::gamma * cc::R * cell.phy.T);
    return production - destruction + gradient_source - compressible;
}

inline void SA::diffusion_SA(cc::face_class& face){
    const double mu = face.otphy.mu;
    const double chi = _chi(face.phy.rho, mu, face.tur.miubl);
    const double mut = face.phy.rho * _fv1(chi) * face.tur.miubl;
    // Reynold 应力
    const double mueff = mut + mu;face.tur.mueff = mueff;
    const double tauxx = mueff*(4.0/3*face.phgrad.ugrad.x - 2.0/3*face.phgrad.vgrad.y);
    const double tauxy = mueff*(face.phgrad.ugrad.y + face.phgrad.vgrad.x);
    const double tauyy = mueff*(4.0/3*face.phgrad.vgrad.y - 2.0/3*face.phgrad.ugrad.x);
    // 热流
    const double lambdaeff = mu/cc::Pr + mut/Prt;
    const cc::vec2 q = -lambdaeff * cc::Cp * face.phgrad.Tgrad;
    // 粘性通量
    const cc::vec4 visF(0,tauxx,tauxy,face.phy.u*tauxx + face.phy.v*tauxy - q.x);
    const cc::vec4 visG(0,tauxy,tauyy,face.phy.u*tauxy + face.phy.v*tauyy - q.y);
    face.visflux.clear();
    face.visflux = face.toface_jacobi(visF,visG);
}

inline void SA::SA_equation_after(cc::cell_class &cell,double dtau,double dt,bool urans){
    double convect = 0.0;double diffusion = 0.0;double new_miubl = 0.0;
    allface(cell){
        cc::face_class *face = cell.faces[i];
        double mu = sutherland::dynamic_viscosity(face->phy.T);
        convect += cell.faces[i]->toface_jacobi(face->phy.u*face->tur.miubl,face->phy.v*face->tur.miubl);
        diffusion += cell.faces[i]->toface_jacobi((mu+face->phy.rho*face->tur.miubl)*face->tur.miublgrad.x
                                                ,(mu+face->phy.rho*face->tur.miubl)*face->tur.miublgrad.x);
    }
    if(!urans) new_miubl = cell.tur.miubl + relax * dtau * ((diffusion - convect)/cell.vol + source_SA(cell));
    else new_miubl = cell.tur.miubl - (1.0/cell.localdt + 1.0/(2*dt))*
                    ((diffusion - convect)/cell.vol + source_SA(cell) + 1/(2*dt)*(cell.tur.miubl-cell.tur.miubl_former));
    cell.tur.miubl = std::min(1e-6,new_miubl);
    // return new_miubl;
}

inline void SA::assemble_visflux(cc::cell_class &cell){
    cell.visflux.clear();
    allface(cell) cell.visflux += cell.faces[i]->visflux * (2*cell.fnorm[i]-1);
}