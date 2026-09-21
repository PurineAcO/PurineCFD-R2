#include "SA.h"
#include "config.h"
#include "physic.h"
#include <cmath>

namespace {

double viscosity_ratio(double rho,double mu,double miubl){
    return rho*(miubl > 0.0 ? miubl : 0.0)/mu;
}

double compute_fv1(double chi){
    return (chi*chi*chi)/(chi*chi*chi + SA::Cv1*SA::Cv1*SA::Cv1);
}

double compute_ft2(double chi){
    return SA::Ct3*std::exp(-SA::Ct4*chi*chi);
}

double compute_fv2(double chi){
    return 1 - chi/(1 + chi*compute_fv1(chi));
}

double compute_g(double r){
    return r + SA::Cw2*(r*r*r*r*r*r - r);
}

double compute_fw(double g){
    constexpr double cw3_squared = SA::Cw3*SA::Cw3;
    constexpr double cw3_sixth = cw3_squared*cw3_squared*cw3_squared;
    const double g_squared = g*g;
    const double g_sixth = g_squared*g_squared*g_squared;
    return g*std::pow((1 + cw3_sixth)/(cw3_sixth + g_sixth),1.0/6);
}

double source_SA(const cc::cell_class& cell){
    double mu = sutherland::dynamic_viscosity(cell.phy.T);
    double chi = viscosity_ratio(cell.phy.rho,mu,cell.tur.miubl);
    double ft2 = compute_ft2(chi);
    double fv2 = compute_fv2(chi);
    double vorticity = std::abs(cell.phy.vgrad.x - cell.phy.ugrad.y);
    const double scaled_nu = cell.tur.miubl*(1.0/(cell.tur.sad*cell.tur.sad))/(SA::kappa*SA::kappa);
    double modified_vorticity =
        std::max(vorticity + fv2*scaled_nu,std::max(0.3*vorticity,1e-20));
    double production = SA::Cb1*(1 - ft2)*modified_vorticity*cell.phy.rho*cell.tur.miubl;
    double r = std::min(scaled_nu/modified_vorticity,SA::rmax);
    double g = compute_g(r);
    double destruction = cell.phy.rho*
                         (SA::Cw1*compute_fw(g) - SA::Cb1/SA::kappa/SA::kappa*ft2)*
                         cell.tur.miubl*cell.tur.miubl*(1.0/(cell.tur.sad*cell.tur.sad));
    double gradient_source = SA::Cb2*SA::inv_sigma*cell.phy.rho*
                             cc::dot(cell.tur.miublgrad,cell.tur.miublgrad);
    // 部分论文中引入了可压缩性修正
    double S2 = 2*cell.phy.ugrad.x*cell.phy.ugrad.x + 2*cell.phy.vgrad.y*cell.phy.vgrad.y + 
                (cell.phy.ugrad.y + cell.phy.vgrad.x)*(cell.phy.ugrad.y + cell.phy.vgrad.x);
    double compressible = SA::C5 * cell.phy.rho * cell.tur.miubl * cell.tur.miubl * S2 / (cc::gamma * cc::R * cell.phy.T);
    return production - destruction + gradient_source - compressible;
}

}

void SA::diffusion_SA(cc::face_class& face){
    const double mu = sutherland::dynamic_viscosity(face.phy.T);
    face.volflux = face.phy.u*face.nor.x + face.phy.v*face.nor.y;
    face.lam = std::abs(face.volflux) + face.phy.a*face.len;
    face.coef = face.tur.miubl + mu/face.phy.rho;
    face.turflux = face.tur.miublgrad.x*face.nor.x + face.tur.miublgrad.y*face.nor.y;
    const double chi = viscosity_ratio(face.phy.rho,mu,face.tur.miubl);
    const double mut = face.phy.rho*compute_fv1(chi)*face.tur.miubl;
    const double mu_eff = mut + mu;
    const double tau_xx = mu_eff*(4.0/3*face.phy.ugrad.x - 2.0/3*face.phy.vgrad.y);
    const double tau_yy = mu_eff*(4.0/3*face.phy.vgrad.y - 2.0/3*face.phy.ugrad.x);
    const double tau_xy = mu_eff*(face.phy.ugrad.y + face.phy.vgrad.x);
    const double lambda_eff = mu/cc::Pr + mut/SA::Prt;
    const cc::vec2 q = {-lambda_eff*cc::Cp*face.phy.Tgrad.x,
                        -lambda_eff*cc::Cp*face.phy.Tgrad.y};
    face.visflux[0] = 0.0;
    face.visflux[1] = tau_xx*face.nor.x + tau_xy*face.nor.y;
    face.visflux[2] = tau_xy*face.nor.x + tau_yy*face.nor.y;
    face.visflux[3] = (face.phy.u*tau_xx + face.phy.v*tau_xy - q.x)*face.nor.x +
                      (face.phy.u*tau_xy + face.phy.v*tau_yy - q.y)*face.nor.y;
}

void SA::SA_equation_RK(cc::cell_class& cell,double rk){
    double rhs = 0.0;
    double velocity_divergence = 0.0;
    for(int i=0;i<cell.ecnt;i++){
        cc::face_class* face = cell.faces[i];
        int outward_sign = 2*cell.fnorm[i] - 1;
        double outward_volume_flux = outward_sign*face->volflux;
        double upwind = cell.tur.miubl;
        if(outward_volume_flux < 0.0){
            if(cell.nei[i] != nullptr){
                upwind = cell.nei[i]->tur.miubl;
            }else{
                upwind = face->tur.miubl;
            }
        }
        double diffusivity = face->coef;
        double outward_gradient_flux = outward_sign*face->turflux;
        rhs += (outward_volume_flux*upwind -
                SA::inv_sigma*diffusivity*outward_gradient_flux)*cell.invvol;
        velocity_divergence += outward_volume_flux*cell.invvol;
    }
    rhs -= source_SA(cell)/cell.phy.rho + cell.tur.miubl*velocity_divergence;
    double next = cell.tur.miubl_former - SA::relax*rk*cell.localdt*rhs;
    cell.tur.miubl_next = std::isfinite(next) ? std::max(next,0.0) : next;
}
