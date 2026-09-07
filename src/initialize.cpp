#include "initialize.h"
#include "classconfig.h"
#include "config.h"
#include "physic.h"
#include <cstdio>

void std_initialize(){
    // 外流: 全场按远场自由流(含攻角)初始化; 无远场配置时回退到速度入口
    const cc::VIL_condition& fs = (cc::FAR_DEFINE.T > 0.0) ? cc::FAR_DEFINE : cc::VIL_DEFINE;
    double rho_inf = fs.p/(cc::R*fs.T);
    double miubl_inf = 3.0 * sutherland::sutherland(fs.T)/rho_inf;
    for(cc::cell_class& cell : cc::CellList){
        cell.phy.T = fs.T;
        cell.phy.p = fs.p;
        cell.phy.u = fs.u;
        cell.phy.v = fs.v;
        cell.phy.rho = cell.phy.p/cell.phy.T/cc::R;
        cell.phy.a = get_sonic_velocity(cell.phy.T);
        cell.phy.e = cc::Cv * cell.phy.T + 0.5*(cell.phy.u*cell.phy.u+cell.phy.v*cell.phy.v);
        cell.tur.miubl = miubl_inf;
    }
    for(cc::face_class& face : cc::FaceList){
        face.tur.miubl = (face.type == cc::WALL) ? 0.0 : miubl_inf;
        if(face.type == cc::INTER){face.face_physic_mid();}
    }
    printf("STD Initialization OK!, u is: %f", cc::CellList[0].phy.u);
}
