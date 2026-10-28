#pragma once

#include "classconfig.hpp"
#include "config.hpp"
#include "physic.hpp"
#include <cstdio>

// 标准初始化
void std_initialize();

inline void std_initialize() {
    // double rho_inf = cc::FAR_DEFINE.phy.rho;
    double miubl_inf = config::miubl_inf_ratio * sutherland::dynamic_viscosity(cc::FAR_DEFINE.phy.T) /
                       cc::FAR_DEFINE.phy.rho;
    allcell {
        icell(i).phy = cc::FAR_DEFINE.phy;
        icell(i).form_otherphy();
        icell(i).tur.miubl = miubl_inf;
    }
    for (cc::face_class& face : cc::FaceList) {
        face.tur.miubl = (face.type == cc::WALL) ? 0.0 : miubl_inf;
        if (face.type == cc::INTER) face.face_physic_mid();
    }
    printf("STD Initialization OK!, u is: %f\n", cc::CellList[0].phy.u);
}
