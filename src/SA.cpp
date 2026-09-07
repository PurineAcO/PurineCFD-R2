#include "SA.h"
#include "classconfig.h"
#include "config.h"
#include "physic.h"
#include <cmath>

// 计算湍流粘度比chi
static double chi_(double rho,double mu,double miubl){return rho * (miubl>0.0?miubl:0.0) / mu;}
// 计算粘度阻尼函数fv1
static double fv1_(double chi){return (chi*chi*chi)/(chi*chi*chi + SA::Cv1*SA::Cv1*SA::Cv1);}
// 计算生产项修正函数ft2
inline double ft2_(double chi){return SA::Ct3 * std::exp(-SA::Ct4 * chi *chi);}
// 计算涡量修正函数fv2
inline double fv2_(double chi){return 1-chi/(1+chi*fv1_(chi));}
// 计算粘性张量S
inline double S_(cc::vec2 ugrad,cc::vec2 vgrad)
{return std::sqrt(2*ugrad.x*ugrad.x + 2*vgrad.y*vgrad.y + (ugrad.y + vgrad.x)*(ugrad.y + vgrad.x));}
// 计算修正涡量Sbl
inline double Sbl_(double Omega,double fv2,double miubl,double sad){return Omega+(fv2*miubl)/(SA::kappa*SA::kappa*sad*sad);}
// 计算无量纲壁面距离r
inline double r_(double miubl,double Sbl,double sad){return miubl/Sbl/SA::kappa/SA::kappa/sad/sad;}
// 计算壁面指标g
inline double g_(double r){return r + SA::Cw2*(r*r*r*r*r*r-r);}
// 计算壁面阻尼函数fw
inline double fw_(double g){return g * pow((1+pow(SA::Cw3,6))/(pow(SA::Cw3,6) + pow(g,6)),1.0/6);}


void SA::diffusion_SA(cc::cell_class &cell){
    for(int j=0;j<4;j++){cell.tur.Ft[j] = 0.0;}
    if(!cc::viscous){ return; }
    for(int i=0;i<cell.ecnt;i++){
        cc::face_class* face = cell.nei[i];
        face->tur.miubl = face->tur.miubl>0 ? face->tur.miubl : 0.0;
        short outer = 2*cell.fnorm[i] - 1;
        double mu = sutherland::sutherland(face->phy.T);
        double chi = chi_(face->phy.rho, mu, face->tur.miubl);
        double fv1 = fv1_(chi);
        double mut = face->phy.rho * fv1 * face->tur.miubl;
        double mu_eff = mut + mu;
        double tau_xx = mu_eff * (4.0/3 * face->phy.ugrad.x - 2.0/3 * face->phy.vgrad.y);
        double tau_yy = mu_eff * (4.0/3 * face->phy.vgrad.y - 2.0/3 * face->phy.ugrad.x);
        double tau_xy = mu_eff * (face->phy.ugrad.y + face->phy.vgrad.x);
        double lambda_eff = mu/cc::Pr + mut/Prt;
        cc::vec2 q ={-lambda_eff*cc::Cp*face->phy.Tgrad.x , -lambda_eff*cc::Cp*face->phy.Tgrad.y};
        cell.tur.Ft[1] += outer * (tau_xx * face->nor.x + tau_xy * face->nor.y);
        cell.tur.Ft[2] += outer * (tau_xy * face->nor.x + tau_yy * face->nor.y);
        cell.tur.Ft[3] += outer * ((face->phy.u*tau_xx + face->phy.v*tau_xy -q.x )*face->nor.x +
                                   (face->phy.u*tau_xy + face->phy.v*tau_yy -q.y )*face->nor.y);
    }
}

double SA::source_SA(cc::cell_class &cell){
    double mu = sutherland::sutherland(cell.phy.T);
    double chi = chi_(cell.phy.rho, mu, cell.tur.miubl);
    double ft2 = ft2_(chi);
    double fv2 = fv2_(chi);
    double Omega = std::abs(cell.phy.vgrad.x - cell.phy.ugrad.y);
    double Sbl = Sbl_(Omega, fv2, cell.tur.miubl, cell.tur.sad);
    double P = Cb1 * (1-ft2) * Sbl * cell.phy.rho * cell.tur.miubl;
    // double S_ = fv3 * S(cell.phy.ugrad,cell.phy.vgrad); // 未启用可压缩时被弃用
    double r = std::min(r_(cell.tur.miubl,Sbl,cell.tur.sad),rmax);
    double g = g_(r);
    double D = cell.phy.rho * (Cw1*fw_(g) - Cb1/kappa/kappa * ft2) * (cell.tur.miubl/cell.tur.sad)*(cell.tur.miubl/cell.tur.sad);
    double G = Cb2 * inv_sigma * cell.phy.rho * cc::dot(cell.tur.miublgrad,cell.tur.miublgrad);
    return P-D+G ;
}

void SA::SA_equation_RK(cc::cell_class &cell, double rk){
    if(!cc::viscous){ return; }
    double res = 0.0;
    for(int j=0;j<cell.ecnt;j++){
        cc::face_class *face = cell.nei[j];
        short outer = 2 * cell.fnorm[j] - 1;
        double mu = sutherland::sutherland(face->phy.T);
        // 一阶迎风对流: mflx(向外为正)>0 上游为本格, <0 入流取邻格或边界给定ν̃
        double mflx = outer * (face->phy.u * face->nor.x + face->phy.v * face->nor.y);
        double nu_up = cell.tur.miubl;
        if(mflx < 0.0){
            if(face->type == cc::INTER){
                cc::cell_class* other = (face->nei[0]->index == cell.index) ? face->nei[1] : face->nei[0];
                nu_up = other->tur.miubl;
            }else{
                nu_up = face->tur.miubl;
            }
        }
        double K = face->phy.rho * face->tur.miubl + mu;
        double gradout = outer * (face->tur.miublgrad.x * face->nor.x + face->tur.miublgrad.y * face->nor.y);
        res += (face->phy.rho * mflx * nu_up - inv_sigma * K * gradout) / cell.vol;
    }
    res -= source_SA(cell);
    double nu = cell.tur.miubl - SA::relax*rk*cell.localdt*res/cell.phy.rho;
    cell.tur.miubl = nu > 0.0 ? nu : 0.0;
}