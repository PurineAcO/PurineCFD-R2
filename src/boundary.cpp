#include "boundary.h"
#include "classconfig.h"
#include "physic.h"
#include "udf.h"
#include <cmath>

void slip_wall_boundary(){
    for(cc::face_class* wall : cc::WallFaces){
        cc::cell_class* c = cc::boundary_findcell(wall);
        wall->phy.u = 0.0;
        wall->phy.v = 0.0;
        wall->phy.T = c->phy.T;
        wall->phy.rho = c->phy.rho;
        wall->tur.miubl = 0.0;
    }
}

void far_field_boundary(){
    double rho_inf = FAR_DEFINE.p/(cc::R*FAR_DEFINE.T);
    double a_inf = get_sonic_velocity(FAR_DEFINE.T);
    double miubl_inf = 3.0*sutherland::dynamic_viscosity(FAR_DEFINE.T)/rho_inf;
    for(cc::face_class* far : cc::FarFaces){
        cc::cell_class* c = cc::boundary_findcell(far);
        double nx = far->nor.x,ny = far->nor.y;
        double len = std::sqrt(nx*nx + ny*ny);
        nx /= len;
        ny /= len;
        if(nx*(far->mid.x - c->center.x) + ny*(far->mid.y - c->center.y) < 0){
            nx = -nx;
            ny = -ny;
        }
        double a = get_sonic_velocity(c->phy.T);
        double vn = c->phy.u*nx + c->phy.v*ny;
        double vt = -1*c->phy.u*ny + c->phy.v*nx;
        double vn_inf = FAR_DEFINE.u*nx + FAR_DEFINE.v*ny;
        double vt_inf = -FAR_DEFINE.u*ny + FAR_DEFINE.v*nx;
        double Rp = vn + 2.0*a/(cc::gamma - 1.0);
        double Rm = vn_inf - 2.0*a_inf/(cc::gamma - 1.0);
        double vn_star = 0.5*(Rp + Rm);
        double a_star = 0.25*(cc::gamma - 1.0)*(Rp - Rm);
        double s,vt_star;
        if(vn_star >= 0.0){
            s = c->phy.p/std::pow(c->phy.rho,cc::gamma);
            vt_star = vt;
        }else{
            s = FAR_DEFINE.p/std::pow(rho_inf,cc::gamma);
            vt_star = vt_inf;
        }
        far->phy.u = vn_star*nx - vt_star*ny;
        far->phy.v = vn_star*ny + vt_star*nx;
        far->phy.rho = std::pow(a_star*a_star/(cc::gamma*s),1.0/(cc::gamma - 1.0));
        far->phy.p = s*std::pow(far->phy.rho,cc::gamma);
        far->phy.T = far->phy.p/(cc::R*far->phy.rho);
        far->tur.miubl = miubl_inf;
    }
}
