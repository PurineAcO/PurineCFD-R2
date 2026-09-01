#include "initialize.h"
#include "classconfig.h"
#include "config.h"
#include "physic.h"
#include "boundary.h"
#include "SA.h"
#include <cstdio>

void std_initialize(){
    cc::total_time = 0.0;
    // 第一步: 初始化所有单元的原始物理量
    for(cc::cell_class& cell : cc::CellList){
        cell.phy.T = cc::VIL_DEFINE.T;
        cell.phy.p = cc::VIL_DEFINE.p;
        cell.phy.u = cc::VIL_DEFINE.u;
        cell.phy.v = cc::VIL_DEFINE.v;
        cell.phy.rho = cell.phy.p/cell.phy.T/cc::R;
        cell.phy.mu = SutherLand(cell.phy.T);
        cell.phy.a = get_sonic_velocity(cell.phy.T);
        cell.tur.miubl = cell.phy.mu/cell.phy.rho * 0.1;
        double chi_ = cell.phy.rho * cell.tur.miubl / cell.phy.mu;
        double chi3 = chi_*chi_*chi_;
        double fv1_ = chi3/(chi3 + SA::Cv1*SA::Cv1*SA::Cv1);
        cell.tur.mut = cell.phy.rho * cell.tur.miubl * fv1_;
        cell.phy.e = get_energy(cell.phy);
    }
    // 第二步: 所有单元初始化完成后, 再统一建构内部面的基本量(避免读到未初始化的邻接单元)
    for(cc::face_class& face : cc::FaceList){
        if(face.type == cc::INTER){face.face_physic_mid();}
    }
    printf("STD Initialization OK!,u is:%f",cc::CellList[0].phy.u);
}