#pragma once

#include "classconfig.hpp"
#include "config.hpp"
#include "physic.hpp"
#include <algorithm>
#include <cmath>

/*
根据Fluent的求解手段,miubl将会在一个内迭代步以后再进行更新,而不是跟随无粘通量进行计算
不论是显式还是隐式,都遵循上述规律.为此,本文件中不仅定义了SA本体的计算还定义了Jacobi的计算
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

// 隐式部分向量
inline std::vector<double> dnu;    // miubl变化量
inline std::vector<double> dnuold; // 上一步miubl变化量
inline std::vector<double> nurhs;  // miubl残差

// 形成粘性通量
void diffusion_SA(cc::face_class& face);
// 组装粘性通量
void assemble_visflux(cc::cell_class& cell);
// 形成湍流中间变量
void build_SA(cc::cell_class& c);
void build_SA(cc::face_class& f);
// 形成源项
double source_SA(cc::cell_class& cell);
// RK子步内求解湍流方程(非守恒形式),目前被弃用
void SA_equation_after(cc::cell_class& cell, double dtau, double dt, bool urans);
// 粘性通量Jacobi矩阵前体B,冻结梯度
std::pair<cc::mat5, cc::mat5> form_pfpgq_diffusion(const cc::face_class& face);
// 将粘性通量Jacobi矩阵前体B转换成B
cc::mat5 form_pfpw_from_pfpgq(const cc::mat5& J,const cc::face_class& face);
// 粘性通量Jacobi矩阵T,冻结非梯度量
cc::mat5 diffusion_grad_jac(const cc::face_class& face,double dn);
// 面上有效扩散系数
double Gamma_eff(cc::cell_class& c,cc::face_class *f);
// 源项Jacobi
void source_jacobi_diag_SA(cc::cell_class& cell,double inv_dt,double cfl);
// 清空向量
void source_vector_clear();
// 分配向量内存
void malloc_vector();
// 形成SA输运方程右侧残差
double right_residual(cc::cell_class& cell);
// SA隐式求解核心
void solvekernel(int c,int lo,int hi);
// 复制SA求解向量
void nucopy();
// 松弛SA求解器
void relaxer();

} // namespace SA

inline void SA::build_SA(cc::cell_class &c){
    c.otphy.mu = sutherland::dynamic_viscosity(c.phy.T);
    c.tur.chi =  c.phy.rho * (c.tur.miubl > 0.0 ? c.tur.miubl : 0.0) / c.otphy.mu;
    c.tur.fv1 = (c.tur.chi*c.tur.chi*c.tur.chi) /(c.tur.chi*c.tur.chi*c.tur.chi + Cv1 * Cv1 * Cv1);
    c.tur.ft2 = Ct3 * std::exp(-SA::Ct4 * c.tur.chi * c.tur.chi);
    c.tur.fv2 =  1 - c.tur.chi / (1 + c.tur.chi * c.tur.fv1);
    c.tur.Omega = std::abs(c.phgrad.vgrad.x - c.phgrad.ugrad.y);
    c.tur.Sbl = std::max(c.tur.Omega + c.tur.fv2 * c.tur.miubl * (1.0 / (c.tur.sad * c.tur.sad)) / (kappa * kappa)
                , std::max(0.3 * c.tur.Omega, 1e-20));
    c.tur.r = std::min(c.tur.miubl * (1.0 / (c.tur.sad * c.tur.sad)) / (kappa * kappa) / c.tur.Sbl, SA::rmax);
    c.tur.g = c.tur.r + SA::Cw2 * (c.tur.r * c.tur.r * c.tur.r * c.tur.r * c.tur.r * c.tur.r - c.tur.r);
    c.tur.fw = c.tur.g * std::pow((1 + pow(Cw3,6)) / (pow(Cw3,6) + pow(c.tur.g, 6)), 1.0 / 6);
    c.tur.S2 = 2 * c.phgrad.ugrad.x * c.phgrad.ugrad.x + 2 * c.phgrad.vgrad.y * c.phgrad.vgrad.y +
        (c.phgrad.ugrad.y + c.phgrad.vgrad.x) * (c.phgrad.ugrad.y + c.phgrad.vgrad.x);
}

inline void SA::build_SA(cc::face_class& f){
    f.otphy.mu = sutherland::dynamic_viscosity(f.phy.T);
    f.tur.chi =  f.phy.rho * (f.tur.miubl > 0.0 ? f.tur.miubl : 0.0) / f.otphy.mu;
    f.tur.fv1 = (f.tur.chi*f.tur.chi*f.tur.chi) /(f.tur.chi*f.tur.chi*f.tur.chi + Cv1 * Cv1 * Cv1);
}

inline void SA::diffusion_SA(cc::face_class& face) {
    build_SA(face);
    // Reynold 应力
    const double mueff = face.phy.rho*face.tur.fv1*face.tur.miubl + face.otphy.mu;
    face.tur.mueff = mueff;
    const double tauxx = mueff * (4.0 / 3 * face.phgrad.ugrad.x - 2.0 / 3 * face.phgrad.vgrad.y);
    const double tauxy = mueff * (face.phgrad.ugrad.y + face.phgrad.vgrad.x);
    const double tauyy = mueff * (4.0 / 3 * face.phgrad.vgrad.y - 2.0 / 3 * face.phgrad.ugrad.x);
    // 热流
    const double lambdaeff = face.otphy.mu / cc::Pr + face.phy.rho*face.tur.fv1*face.tur.miubl / Prt;
    const cc::vec2 q = -lambdaeff * cc::Cp * face.phgrad.Tgrad;
    // 粘性通量
    const cc::vec4 visF(0, tauxx, tauxy, face.phy.u * tauxx + face.phy.v * tauxy - q.x);
    const cc::vec4 visG(0, tauxy, tauyy, face.phy.u * tauxy + face.phy.v * tauyy - q.y);
    face.visflux.clear();
    face.visflux = face.toface_jacobi(visF, visG);
}

inline double SA::source_SA(cc::cell_class& cell) {
    build_SA(cell);
    double production = Cb1 * (1 - cell.tur.ft2) * cell.tur.Sbl * cell.phy.rho * cell.tur.miubl;
    double destruction = cell.phy.rho * (Cw1 * cell.tur.fw - Cb1 / kappa / kappa * cell.tur.ft2) *
                         cell.tur.miubl * cell.tur.miubl * (1.0 / (cell.tur.sad * cell.tur.sad));
    // 部分论文中引入了可压缩性修正
    double compressible = C5 * cell.phy.rho * cell.tur.miubl * cell.tur.miubl * cell.tur.S2 /
                          (cc::gamma * cc::R * cell.phy.T);
    return production - destruction - compressible;
}


inline std::pair<cc::mat5, cc::mat5> SA::form_pfpgq_diffusion(const cc::face_class& face){
    const cc::vecgrad& g = face.phgrad;
    const cc::vec2& bg = face.tur.miublgrad;
    const double rho = face.phy.rho, u = face.phy.u, v = face.phy.v;
    const double fv1 = face.tur.fv1;
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

inline double SA::Gamma_eff(cc::cell_class& c, cc::face_class *f){
    if(f->type == cc::WALL) return inv_sigma * f->otphy.mu;
    else if(f->type == cc::FAR){return inv_sigma*(f->otphy.mu + f->phy.rho *f->tur.miubl);}
    else if(f->type == cc::INTER){
        return (1+Cb2)*inv_sigma * (f->otphy.mu + f->phy.rho * f->tur.miubl) - 
            Cb2*inv_sigma * (c.otphy.mu + c.phy.rho*c.tur.miubl);}
    else return 0.0;
}

// inline void SA::source_vector_clear(){
// #pragma omp parallel for schedule(static)
//     for(int i=0;i<(int) dnu.size();i++) dnu[i]=0.0;
// }

inline void SA::source_jacobi_diag_SA(cc::cell_class &cell,double inv_dt,double cfl){
    cell.tur.diag = cell.vol * cell.phy.rho * inv_dt;
    cell.tur.diag += cell.phy.rho * cell.lambda / cfl;
    allface(cell){
        cc::face_class* f = cell.faces[i];
        const double mdot = (f->nei[0] == &cell ? 1.0 : -1.0) * f->convect.c;
        cell.tur.diag += std::max(mdot, 0.0);
        cell.tur.diag += Gamma_eff(cell,f) * f->len / normal_dist(f);
    }
    double chidft2 = -2*Ct4 *cell.tur.chi*cell.tur.chi *cell.tur.ft2;
    double chidfv2 = -cell.tur.chi * (1-3*cell.tur.chi*cell.tur.fv1*(1-cell.tur.fv1))/pow(1+cell.tur.chi*cell.tur.fv1,2);
    double cw3_6 = pow(Cw3,6);
    double g6 = pow(cell.tur.g,6);
    double h = pow((1+cw3_6)/(cw3_6+g6),1.0/6);
    double rdfw = cell.tur.r * h * cw3_6/(cw3_6+g6) *
                   (1 + Cw2*(6*pow(cell.tur.r,5) - 1));
    double rhoS = cell.phy.rho * cell.tur.Sbl;
    cell.tur.diag -= cell.vol * ( rhoS * Cb1 * (1-cell.tur.ft2-chidft2)
                    + rhoS * Cb1 * cell.tur.r * ((cell.tur.fv2+chidfv2)*(1-cell.tur.ft2) +2*cell.tur.ft2 + chidft2 )
                    - rhoS * Cw1 * kappa*kappa*cell.tur.r *(2*cell.tur.fw+rdfw-cell.tur.r*rdfw*(cell.tur.fv2+chidfv2)) );
    cell.tur.diag += cell.vol * 2*C5*cell.phy.rho*cell.tur.miubl*cell.tur.S2/(cc::gamma*cc::R*cell.phy.T);
}

inline void SA::malloc_vector(){
    if((int) dnu.size() == cc::cell_num)return;
    dnu.resize(cc::cell_num);dnuold.resize(cc::cell_num);nurhs.resize(cc::cell_num);
}


inline double SA::right_residual(cc::cell_class& cell){
    double res = 0.0;
    allface(cell){
        cc::face_class* f = cell.faces[i];
        const cc::cell_class* nb = cc::other_cell(f,&cell);
        const double nu_o = nb ? nb->tur.miubl : f->tur.miubl;
        const double mdot = (f->nei[0] == &cell ? 1.0 : -1.0) * f->convect.c;
        res += mdot * (mdot > 0.0 ? cell.tur.miubl : nu_o);
        res += Gamma_eff(cell, f)*f->len/normal_dist(f)*(cell.tur.miubl - nu_o);
    }
    res -= cell.vol * source_SA(cell);
    return res;
}

inline void SA::solvekernel(int c, int lo, int hi){
    cc::cell_class& cell = icell(c);
    double acc = nurhs[c];
    allface(cell){
        cc::face_class* f = cell.faces[i];
        if(f->type != cc::INTER) continue;
        cc::cell_class* nb = cc::other_cell(f, &cell);
        const int j = nb->index-1;
        const double dj = (j>=lo && j<hi) ? dnu[j] : dnuold[j];
        const double mdot = (f->nei[0] == &cell ? 1.0 : -1.0) * f->convect.c;
        acc += (std::min(mdot, 0.0) - Gamma_eff(cell, f)*f->len/normal_dist(f))*dj;
    }
    dnu[c] = -acc/cell.tur.diag;
}

inline void SA::nucopy(){
#pragma omp parallel for schedule(static)
    for(int i=0;i<(int) dnu.size();i++) dnuold[i] = dnu[i];
}

inline void SA::relaxer(){
#pragma omp parallel for schedule(static)
    allcell{
        double delta_nu = dnu[i];
        if(!std::isfinite(dnu[i])) {delta_nu = 0.0;}
        icell(i).tur.miubl = std::max(0.0,icell(i).tur.miubl + delta_nu);
    }
}