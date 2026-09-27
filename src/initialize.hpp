#pragma once

#include "classconfig.hpp"
#include "config.hpp"
#include "physic.hpp"
#include "udf.hpp"
#include <cstdio>

// 标准初始化
void std_initialize();

inline void std_initialize(){
    double rho_inf = FAR_DEFINE.p/(cc::R*FAR_DEFINE.T);
    double miubl_inf = 3.0*sutherland::dynamic_viscosity(FAR_DEFINE.T)/rho_inf;
    for(cc::cell_class& cell : cc::CellList){
        cell.phy = cc::vecp(rho_inf,FAR_DEFINE.u,FAR_DEFINE.v,FAR_DEFINE.T);
        cell.otphy = cc::otphy(get_sonic_velocity(FAR_DEFINE.T),FAR_DEFINE.p,get_energy(cell.phy));
        cell.tur.miubl = miubl_inf;
    }
    for(cc::face_class& face : cc::FaceList){
        face.tur.miubl = (face.type == cc::WALL) ? 0.0 : miubl_inf;
        if(face.type == cc::INTER) face.face_physic_mid();
    }
    printf("STD Initialization OK!, u is: %f\n",cc::CellList[0].phy.u);
}
