#include "boundary.h"
#include "classconfig.h"
#include "config.h"
#include "physic.h"
#include <cmath>

// 有滑移不穿透壁面
void slip_wall_boundary(){
    for(cc::face_class* wall : cc::WallFaces){

        cc::cell_class* c = cc::boundary_findcell(wall);
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
    double a_inf = get_sonic_velocity(cc::FAR_DEFINE.T);
    for(cc::face_class* far : cc::FARFaces){

        cc::cell_class* c = cc::boundary_findcell(far);
        double nx = far->nor.x, ny = far->nor.y;double len = std::sqrt(nx*nx+ny*ny);nx /= len; ny /= len;
        if(nx*(far->mid.x - c->center.x) + ny*(far->mid.y - c->center.y) < 0){ nx = -nx; ny = -ny; }

        // Riemann 不变量计算
        double a = get_sonic_velocity(c->phy.T);
        double vn = c->phy.u*nx + c->phy.v*ny;
        double vt = -1*c->phy.u*ny + c->phy.v*nx;
        double vn_inf = cc::FAR_DEFINE.u*nx + cc::FAR_DEFINE.v*ny;
        double vt_inf = -cc::FAR_DEFINE.u*ny + cc::FAR_DEFINE.v*nx;
        double Rp = vn + 2.0*a/(cc::gamma-1.0);                    // 第一不变量,出波
        double Rm = vn_inf - 2.0*a_inf/(cc::gamma-1.0);            // 第二不变量,入波
        double vn_star = 0.5*(Rp + Rm);
        double a_star = 0.25*(cc::gamma-1.0)*(Rp - Rm);
        double s, vt_star;
        if(vn_star >= 0.0){ // 出流,取内部
            s = c->phy.p/std::pow(c->phy.rho, cc::gamma);  // 第三不变量,熵的衍生物
            vt_star = vt;                                          // 第四不变量,切向速度
        }else{ // 入流,取来流
            s = cc::FAR_DEFINE.p/std::pow(rho_inf, cc::gamma);
            vt_star = vt_inf;
        }

        // 还原物理量
        far->phy.u = vn_star*nx - vt_star*ny;
        far->phy.v = vn_star*ny + vt_star*nx;
        far->phy.rho = std::pow(a_star*a_star/(cc::gamma*s), 1.0/(cc::gamma-1.0));
        far->phy.p = s*std::pow(far->phy.rho, cc::gamma);
        far->phy.T = far->phy.p/(cc::R*far->phy.rho);
    }
}

// 压力远场赋值函数
void farfield_boundary_init(double p, double Ma, double T, double AOA){
    cc::FAR_DEFINE.p = p;
    cc::FAR_DEFINE.T = T;
    cc::FAR_DEFINE.u = cos(deg2rad(AOA)) * Ma * get_sonic_velocity(T);
    cc::FAR_DEFINE.v = sin(deg2rad(AOA)) * Ma * get_sonic_velocity(T);
}
