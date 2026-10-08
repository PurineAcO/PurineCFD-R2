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
double source_SA(const cc::cell_class& cell);
// RK子步内求解湍流方程(非守恒形式)
void SA_equation_after(cc::cell_class& cell, double dtau, double dt, bool urans);
// 湍流粘度比
double _chi(double rho, double mu, double miubl);
// 粘度阻尼函数
double _fv1(double chi);
// 转捩控制函数
double _ft2(double chi);
// 涡量辅助函数
double _fv2(double chi);
// 无量纲距离II
double _g(double r);
// 破坏项控制函数
double _fw(double g);
// 粘性通量Jacobi矩阵前体B,冻结梯度
std::pair<cc::mat5, cc::mat5> form_pfpgq_diffusion(const cc::face_class& face);
// 将粘性通量Jacobi矩阵前体B转换成B
cc::mat5 form_pfpw_from_pfpgq(const cc::mat5& J,const cc::face_class& face);
// 粘性通量Jacobi矩阵T,冻结非梯度量
cc::mat5 diffusion_grad_jac(const cc::face_class& face,double dn);

} // namespace SA

namespace sutherland {

// sutherland粘度对温度的导数
double _MUT(double mu,double T);

} // namespace sutherland



inline double SA::_chi(double rho, double mu, double miubl) {
    return rho * (miubl > 0.0 ? miubl : 0.0) / mu;
}

inline double SA::_fv1(double chi) {
    return (chi * chi * chi) / (chi * chi * chi + SA::Cv1 * SA::Cv1 * SA::Cv1);
}

inline double SA::_ft2(double chi) {
    return SA::Ct3 * std::exp(-SA::Ct4 * chi * chi);
}

inline double SA::_fv2(double chi) {
    return 1 - chi / (1 + chi * SA::_fv1(chi));
}

inline double SA::_g(double r) {
    return r + SA::Cw2 * (r * r * r * r * r * r - r);
}

inline double SA::_fw(double g) {
    constexpr double cw3_squared = SA::Cw3 * SA::Cw3;
    constexpr double cw3_sixth = cw3_squared * cw3_squared * cw3_squared;
    const double g_squared = g * g;
    const double g_sixth = g_squared * g_squared * g_squared;
    return g * std::pow((1 + cw3_sixth) / (cw3_sixth + g_sixth), 1.0 / 6);
}

inline double SA::source_SA(const cc::cell_class& cell) {
    double mu = cell.otphy.mu;
    double chi = _chi(cell.phy.rho, mu, cell.tur.miubl);
    double ft2 = _ft2(chi);
    double fv2 = _fv2(chi);
    double vorticity = std::abs(cell.phgrad.vgrad.x - cell.phgrad.ugrad.y);
    const double scaled_nu =
        cell.tur.miubl * (1.0 / (cell.tur.sad * cell.tur.sad)) / (SA::kappa * SA::kappa);
    double modified_vorticity =
        std::max(vorticity + fv2 * scaled_nu, std::max(0.3 * vorticity, 1e-20));
    double production = SA::Cb1 * (1 - ft2) * modified_vorticity * cell.phy.rho * cell.tur.miubl;
    double r = std::min(scaled_nu / modified_vorticity, SA::rmax);
    double g = _g(r);
    double destruction = cell.phy.rho * (SA::Cw1 * _fw(g) - SA::Cb1 / SA::kappa / SA::kappa * ft2) *
                         cell.tur.miubl * cell.tur.miubl * (1.0 / (cell.tur.sad * cell.tur.sad));
    double gradient_source =
        SA::Cb2 * SA::inv_sigma * cell.phy.rho * cc::dot(cell.tur.miublgrad, cell.tur.miublgrad);
    // 部分论文中引入了可压缩性修正
    double S2 =
        2 * cell.phgrad.ugrad.x * cell.phgrad.ugrad.x +
        2 * cell.phgrad.vgrad.y * cell.phgrad.vgrad.y +
        (cell.phgrad.ugrad.y + cell.phgrad.vgrad.x) * (cell.phgrad.ugrad.y + cell.phgrad.vgrad.x);
    double compressible = SA::C5 * cell.phy.rho * cell.tur.miubl * cell.tur.miubl * S2 /
                          (cc::gamma * cc::R * cell.phy.T);
    return production - destruction + gradient_source - compressible;
}

inline double sutherland::_MUT(double mu,double T){
    return mu*(T+3*sutherland::Ts)/(2*T*(T+sutherland::Ts));
}

inline std::pair<cc::mat5, cc::mat5> SA::form_pfpgq_diffusion(const cc::face_class& face){
    const cc::vecgrad& g = face.phgrad;
    const cc::vec2& bg = face.tur.miublgrad;
    const double rho = face.phy.rho, u = face.phy.u, v = face.phy.v;
    const double fv1 = SA::_fv1(SA::_chi(rho, face.otphy.mu, face.tur.miubl));
    const double D = fv1 + 3*fv1*(1-fv1);
    const double A2 = g.ugrad.x*4/3 - g.vgrad.y*2/3;
    const double A3 = g.ugrad.y + g.vgrad.x;
    const double A4 = g.vgrad.y*4/3 - g.ugrad.x*2/3;
    const double B2 = A2*D, B3 = A3*D, B4 = A4*D;
    const double mueff = face.otphy.mu + rho*fv1*face.tur.miubl;
    const double tauxx = mueff*A2, tauxy = mueff*A3, tauyy = mueff*A4;
    const double MUT = sutherland::_MUT(face.otphy.mu,face.phy.T);
    const double E1 = bg.x*SA::inv_sigma, E2 = bg.y*SA::inv_sigma;
    cc::vec5 c1(0,0,0,0,0);
    cc::vec5 x1(B2*face.tur.miubl,0,0,MUT*A2,B2*rho);
    cc::vec5 y1(B3*face.tur.miubl,0,0,MUT*A3,B3*rho);
    cc::vec5 e1(face.tur.miubl*(u*B2+v*B3+D*cc::Cp/SA::Prt*g.Tgrad.x),
                tauxx,tauxy,MUT*(cc::Cp/cc::Pr*g.Tgrad.x+u*A2+v*A3),
                rho*(u*B2+v*B3+D*cc::Cp/SA::Prt*g.Tgrad.x));
    cc::vec5 t1(face.tur.miubl*E1,0,0,MUT*E1,rho*E1);
    cc::mat5 F = cc::mat5(c1,x1,y1,e1,t1);
    cc::vec5 c2(0,0,0,0,0);
    cc::vec5 x2(B3*face.tur.miubl,0,0,MUT*A3,B3*rho);
    cc::vec5 y2(B4*face.tur.miubl,0,0,MUT*A4,B4*rho);
    cc::vec5 e2(face.tur.miubl*(u*B3+v*B4+D*cc::Cp/SA::Prt*g.Tgrad.y),
                tauxy,tauyy,MUT*(cc::Cp/cc::Pr*g.Tgrad.y+u*A3+v*A4),
                rho*(u*B3+v*B4+D*cc::Cp/SA::Prt*g.Tgrad.y));
    cc::vec5 t2(face.tur.miubl*E2,0,0,MUT*E2,rho*E2);
    cc::mat5 G = cc::mat5(c2,x2,y2,e2,t2);
    return std::pair<cc::mat5,cc::mat5>(F,G);
}

inline cc::mat5 SA::form_pfpw_from_pfpgq(const cc::mat5& J,const cc::face_class& face){
    const double rho = face.phy.rho, u = face.phy.u, v = face.phy.v, T = face.phy.T, nu = face.tur.miubl;
    const double Cv = cc::Cv, ke = (u*u+v*v)/(2*Cv*rho);
    const cc::vec5 n_rho(1,0,0,0,0);
    const cc::vec5 n_u(-u/rho,1/rho,0,0,0);
    const cc::vec5 n_v(-v/rho,0,1/rho,0,0);
    const cc::vec5 n_T(ke-T/rho,-u/(Cv*rho),-v/(Cv*rho),1/(Cv*rho),0);
    const cc::vec5 n_nu(-nu/rho,0,0,0,1/rho);
    cc::mat5 M;
    for(int r=0;r<5;r++){
        const cc::vec5& a = cc::mat5_row(J,r);
        cc::mat5_row(M,r) = a.rho*n_rho+a.u*n_u+a.v*n_v+a.T*n_T+a.miubl*n_nu;
    }
    return M;
}

inline cc::mat5 SA::diffusion_grad_jac(const cc::face_class& face,double dn){
    const double len = face.len, nx = face.nor.x/len, ny = face.nor.y/len;
    const double rho = face.phy.rho, u = face.phy.u, v = face.phy.v, nu = face.tur.miubl, T = face.phy.T;
    const double mu = face.otphy.mu, mut = std::max(face.tur.mueff-mu,0.0), m = face.tur.mueff;
    const double k = cc::Cp*(mu/cc::Pr+mut/SA::Prt), sg = SA::inv_sigma*(mu+rho*nu), Cv = cc::Cv;
    const double Nu[5] = {-u/rho,1/rho,0,0,0};
    const double Nv[5] = {-v/rho,0,1/rho,0,0};
    const double NT[5] = {(u*u+v*v)/(2*Cv*rho)-T/rho,-u/(Cv*rho),-v/(Cv*rho),1/(Cv*rho),0};
    const double Nb[5] = {-nu/rho,0,0,0,1/rho};
    double a[5][5] = {};
    for(int i=0;i<5;i++){
        const double a2 = (4.0/3*nx*Nu[i]-2.0/3*ny*Nv[i])/dn;
        const double a3 = (ny*Nu[i]+nx*Nv[i])/dn;
        const double a4 = (4.0/3*ny*Nv[i]-2.0/3*nx*Nu[i])/dn;
        const double tx = m*(nx*a2+ny*a3), ty = m*(nx*a3+ny*a4);
        a[1][i] = len*tx;
        a[2][i] = len*ty;
        a[3][i] = len*(u*tx+v*ty+k*NT[i]/dn);
        a[4][i] = len*sg*Nb[i]/dn;
    }
    return cc::mat5_of(a);
}

inline void SA::diffusion_SA(cc::face_class& face) {
    const double mu = face.otphy.mu;
    const double chi = _chi(face.phy.rho, mu, face.tur.miubl);
    const double mut = face.phy.rho * _fv1(chi) * face.tur.miubl;
    // Reynold 应力
    const double mueff = mut + mu;
    face.tur.mueff = mueff;
    const double tauxx = mueff * (4.0 / 3 * face.phgrad.ugrad.x - 2.0 / 3 * face.phgrad.vgrad.y);
    const double tauxy = mueff * (face.phgrad.ugrad.y + face.phgrad.vgrad.x);
    const double tauyy = mueff * (4.0 / 3 * face.phgrad.vgrad.y - 2.0 / 3 * face.phgrad.ugrad.x);
    // 热流
    const double lambdaeff = mu / cc::Pr + mut / Prt;
    const cc::vec2 q = -lambdaeff * cc::Cp * face.phgrad.Tgrad;
    // 粘性通量
    const cc::vec4 visF(0, tauxx, tauxy, face.phy.u * tauxx + face.phy.v * tauxy - q.x);
    const cc::vec4 visG(0, tauxy, tauyy, face.phy.u * tauxy + face.phy.v * tauyy - q.y);
    face.visflux.clear();
    face.visflux = face.toface_jacobi(visF, visG);
}

inline void SA::SA_equation_after(cc::cell_class& cell, double dtau, double dt, bool urans) {
    double convect = 0.0;
    double diffusion = 0.0;
    double new_miubl = 0.0;
    allface(cell) {
        cc::face_class* face = cell.faces[i];
        const double mu = sutherland::dynamic_viscosity(face->phy.T);
        // 统一成外法向
        const double outer = 2 * cell.fnorm[i] - 1;
        const double coef = inv_sigma * (mu + face->phy.rho * face->tur.miubl);
        // 对流项必须带 rho, 其余各项已带 rho, 统一除以 rho 才是非守恒形式
        convect += outer * face->toface_jacobi(face->phy.rho * face->phy.u * face->tur.miubl,
                                               face->phy.rho * face->phy.v * face->tur.miubl);
        diffusion +=
            outer * face->toface_jacobi(coef * face->tur.miublgrad.x, coef * face->tur.miublgrad.y);
    }
    const double rhs = (diffusion - convect) / (cell.vol) + source_SA(cell);
    if (!urans) {
        new_miubl = (cell.conserformer.c * cell.tur.miubl + relax * dtau * rhs) / cell.phy.rho;
    } else {
        new_miubl = (cell.conserformer.c * cell.tur.miubl +
                     (rhs - (cell.tur.miubl - cell.tur.miubl_former) / (2 * dt)) /
                         (1.0 / dtau + 1.0 / (2 * dt))) /
                    cell.phy.rho;
    }
    cell.tur.miubl = std::max(0.0, new_miubl);
}

inline void SA::assemble_visflux(cc::cell_class& cell) {
    cell.visflux.clear();
    allface(cell) cell.visflux += cell.faces[i]->visflux * (2 * cell.fnorm[i] - 1);
}