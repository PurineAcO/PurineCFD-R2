#include "initialize.h"
#include "classconfig.h"
#include "config.h"
#include "physic.h"
#include <cstdio>

void std_initialize(){
    cc::total_time = 0.0;
    for(cc::cell_class& cell : cc::CellList){
        cell.phy.T = cc::VIL_DEFINE.T;
        cell.phy.p = cc::VIL_DEFINE.p;
        cell.phy.u = cc::VIL_DEFINE.u;
        cell.phy.v = cc::VIL_DEFINE.v;
        cell.phy.rho = cell.phy.p/cell.phy.T/cc::R;
        cell.phy.a = get_sonic_velocity(cell.phy.T);
        cell.phy.e = cc::Cv * cell.phy.T + 0.5*(cell.phy.u*cell.phy.u+cell.phy.v*cell.phy.v);
    }
    for(cc::face_class& face : cc::FaceList){
        if(face.type == cc::INTER){face.face_physic_mid();}
    }
    printf("STD Initialization OK!, u is: %f", cc::CellList[0].phy.u);
}
