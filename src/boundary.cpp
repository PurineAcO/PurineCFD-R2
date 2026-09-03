#include "classconfig.h"
#include "config.h"
#include <cmath>

void slip_wall_boundary(){
    for(cc::face_class* wall : cc::WallFaces){
        cc::cell_class* c = wall->nei[0] ? wall->nei[0] : wall->nei[1];
        double nx = wall->nor.x, ny = wall->nor.y;
        if(nx*(wall->mid.x - c->center.x) + ny*(wall->mid.y - c->center.y) < 0){ nx = -nx; ny = -ny; }
        double vn = (c->phy.u*nx + c->phy.v*ny) / (nx*nx + ny*ny);
        wall->phy.u = c->phy.u - vn*nx;
        wall->phy.v = c->phy.v - vn*ny;
        wall->phy.T = c->phy.T;
        wall->phy.rho = c->phy.rho;
    }
}

// 压力远场
void far_field_boundary(){
    double rho_inf = cc::FAR_DEFINE.p/(cc::R*cc::FAR_DEFINE.T);
    double a_inf = std::sqrt(cc::gamma*cc::R*cc::FAR_DEFINE.T);
    for(cc::face_class* face : cc::FARFaces){
        cc::cell_class* c = face->nei[0] ? face->nei[0] : face->nei[1];
        if(!c){ continue; }
        double nx = face->nor.x, ny = face->nor.y;
        double len = std::hypot(nx, ny);
        if(len < 1e-30){ continue; }
        nx /= len; ny /= len;
        if(nx*(face->mid.x - c->center.x) + ny*(face->mid.y - c->center.y) < 0){ nx = -nx; ny = -ny; }
        double rho = c->phy.rho, p = c->phy.p, T = c->phy.T;
        double u = c->phy.u, v = c->phy.v;
        double a = std::sqrt(cc::gamma*cc::R*T);
        double vn = u*nx + v*ny;
        double vt = -u*ny + v*nx;
        double vn_inf = cc::FAR_DEFINE.u*nx + cc::FAR_DEFINE.v*ny;
        double vt_inf = -cc::FAR_DEFINE.u*ny + cc::FAR_DEFINE.v*nx;
        double Rp = vn + 2.0*a/(cc::gamma-1.0);
        double Rm = vn_inf - 2.0*a_inf/(cc::gamma-1.0);
        double vn_star = 0.5*(Rp + Rm);
        double a_star = 0.25*(cc::gamma-1.0)*(Rp - Rm);
        double s, vt_star;
        if(vn_star >= 0.0){
            s = p/std::pow(rho, cc::gamma);
            vt_star = vt;
        }else{
            s = cc::FAR_DEFINE.p/std::pow(rho_inf, cc::gamma);
            vt_star = vt_inf;
        }
        double rho_star = std::pow(a_star*a_star/(cc::gamma*s), 1.0/(cc::gamma-1.0));
        double p_star = s*std::pow(rho_star, cc::gamma);
        double u_star = vn_star*nx - vt_star*ny;
        double v_star = vn_star*ny + vt_star*nx;
        face->phy.u = u_star;
        face->phy.v = v_star;
        face->phy.rho = rho_star;
        face->phy.p = p_star;
        face->phy.T = p_star/(cc::R*rho_star);
    }
}
