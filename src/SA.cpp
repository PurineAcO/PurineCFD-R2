#include "classconfig.h"
#include "config.h"
#include "SA.h"
#include <vector>

// 计算湍流粘度比chi
double chi(double miubl, double miu,double rho){return (rho*miubl)/miu; }

// 计算粘度阻尼函数fv1
double fv1(double chi){return (chi*chi*chi)/(chi*chi*chi+SA::Cv1*SA::Cv1*SA::Cv1);}

// S-A扩散项函数
void SA::SA_diffusion(cc::cell_class &cell){
    cell.diffusion = {0.0,0.0,0.0,0.0,0.0};
    for(int i=0;i<cell.ecnt;i++){
        cc::face_class* face = cell.nei[i];bool outer = cell.fnorm[i];
        double chi_ = chi(face->tur.miubl,face->phy.mu,face->phy.rho);
        double fv1_ = fv1(chi_);
        double mu_t = face->phy.rho * face->tur.miubl * fv1_;
        double mu_eff = mu_t + face->phy.mu;
        double lambda_eff = mu_t/SA::Prt + face->phy.mu/SA::Pr;
        double tau_xx = mu_eff * (4.0/3.0 * face->phy.ugrad.x - 2.0/3.0 * face->phy.vgrad.y);
        double tau_xy = mu_eff * (face->phy.ugrad.y + face->phy.vgrad.x);
        double tau_yy = mu_eff * (4.0/3.0 * face->phy.vgrad.y - 2.0/3.0 * face->phy.ugrad.x);
        cc::vec2 q = (lambda_eff * cc::Cp) * face->phy.Tgrad;
        std::vector<double> F = {0,tau_xx,tau_xy,
                                face->phy.u * tau_xx + face->phy.v*tau_xy + q.x,
                                SA::inv_sigma * (face->phy.mu + face->phy.rho * face->tur.miubl) * face->tur.miublgrad.x};
        std::vector<double> G = {0,tau_xy,tau_yy,
                                face->phy.u * tau_xy + face->phy.v * tau_yy + q.y,
                                SA::inv_sigma * (face->phy.mu + face->phy.rho * face->tur.miubl) * face->tur.miublgrad.y};
        for(int i=0;i<5;i++){cell.diffusion[i] += (2*outer-1) * (F[i] * face->nor.x + G[i] * face->nor.y);}
    }
}