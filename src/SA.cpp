#include "classconfig.h"
#include "config.h"
#include "SA.h"
#include <cmath>
#include <cstdlib>
#include <vector>

namespace SA {

// 计算湍流粘度比chi
inline double chi(double miubl, double miu,double rho){return (rho*miubl)/miu; }

// 计算粘度阻尼函数fv1
inline double fv1(double chi){return (chi*chi*chi)/(chi*chi*chi+Cv1*Cv1*Cv1);}

// 计算生产项修正函数ft2
inline double ft2(double chi){return Ct3 * std::exp(-Ct4 * chi *chi);}

// 计算涡量修正函数fv2
inline double fv2(double chi){return 1-chi/(1+chi*fv1(chi));}

// 计算粘性张量S
inline double S(cc::vec2 ugrad,cc::vec2 vgrad)
{return std::sqrt(2*ugrad.x*ugrad.x + 2*vgrad.y*vgrad.y + (ugrad.y + vgrad.x)*(ugrad.y + vgrad.x));}

// 计算修正涡量Sbl
inline double Sbl(double Omega,double fv2,double miubl,double sad){return Omega+(fv2*miubl)/(kappa*kappa*sad*sad);}

// 计算无量纲壁面距离r
inline double r(double miubl,double Sbl,double sad){return miubl/Sbl/kappa/kappa/sad/sad;}

// 计算壁面指标g
inline double g(double r){return r + Cw2*(r*r*r*r*r*r-r);}

// 计算壁面阻尼函数fw
inline double fw(double g){return g * pow((1+pow(Cw3,6))/(pow(Cw3,6) + pow(g,6)),1.0/6);}

// S-A扩散项函数
void SA_diffusion(cc::cell_class &cell){
    cell.diffusion = {0.0,0.0,0.0,0.0,0.0};
    cell.tur.mut = cell.phy.rho * cell.tur.miubl * fv1(chi(cell.tur.miubl,cell.phy.mu,cell.phy.rho));
    for(int i=0;i<cell.ecnt;i++){
        cc::face_class* face = cell.nei[i];bool outer = cell.fnorm[i];
        double chi_ = chi(face->tur.miubl,face->phy.mu,face->phy.rho);
        double fv1_ = fv1(chi_);
        face->tur.mut = face->phy.rho * face->tur.miubl * fv1_;
        double mu_eff = face->tur.mut + face->phy.mu;
        double lambda_eff = face->tur.mut/Prt + face->phy.mu/Pr;
        double tau_xx = mu_eff * (4.0/3.0 * face->phy.ugrad.x - 2.0/3.0 * face->phy.vgrad.y);
        double tau_xy = mu_eff * (face->phy.ugrad.y + face->phy.vgrad.x);
        double tau_yy = mu_eff * (4.0/3.0 * face->phy.vgrad.y - 2.0/3.0 * face->phy.ugrad.x);
        cc::vec2 q = (lambda_eff * cc::Cp) * face->phy.Tgrad;
        std::vector<double> F = {0,tau_xx,tau_xy,
                                face->phy.u * tau_xx + face->phy.v*tau_xy + q.x,
                                inv_sigma * (face->phy.mu + face->phy.rho * face->tur.miubl) * face->tur.miublgrad.x};
        std::vector<double> G = {0,tau_xy,tau_yy,
                                face->phy.u * tau_xy + face->phy.v * tau_yy + q.y,
                                inv_sigma * (face->phy.mu + face->phy.rho * face->tur.miubl) * face->tur.miublgrad.y};
        for(int j=0;j<5;j++){cell.diffusion[j] += (2*outer-1) * (F[j] * face->nor.x + G[j] * face->nor.y);}
    }
}

// S-A源项函数
void SA_source(cc::cell_class &cell){
    cell.source = {0.0,0.0,0.0,0.0,0.0};
    double ft2_ = ft2(chi(cell.tur.miubl,cell.phy.mu,cell.phy.rho));
    double fv2_ = fv2(chi(cell.tur.miubl,cell.phy.mu,cell.phy.rho));
    double Omega = std::abs(cell.phy.vgrad.x - cell.phy.ugrad.y);
    double Sbl_ = Sbl(Omega, fv2_, cell.tur.miubl, cell.sad);
    double P = Cb1 * (1-ft2_) * Sbl_ * cell.phy.rho * cell.tur.miubl;
    // double S_ = fv3 * S(cell.phy.ugrad,cell.phy.vgrad); // 未启用可压缩时被弃用
    double r_ = std::min(r(cell.tur.miubl,Sbl_,cell.sad),rmax);
    double g_ = g(r_);
    double D = cell.phy.rho * (Cw1*fw(g_) - Cb1/kappa/kappa * ft2_) * (cell.tur.miubl/cell.sad)*(cell.tur.miubl/cell.sad);
    double G = Cb2 * inv_sigma * cell.phy.rho * cc::dot(cell.tur.miublgrad,cell.tur.miublgrad);
    cell.source = {0.0,0.0,0.0,0.0,P-D+G};
}

}