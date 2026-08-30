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
        cell.phy.mu = SutherLand(cell.phy.T);
        cell.phy.a = get_sonic_velocity(cell.phy.T);
        cell.tur.miubl = cell.phy.mu/cell.phy.rho * 0.1;
        cell.phy.e = get_energy(cell.phy);
    }
    printf("STD Initialization OK!,u is:%f",cc::CellList[0].phy.u);
}