#include "classconfig.h"
#include "config.h"
#include "timarch.h"
#include "SA.h"
#include <cstdlib>
#include <vector>

void local_timestep(cc::cell_class &cell){
    std::vector<double> ds = {0.0,0.0};
    for(int i=0;i<cell.ecnt;i++){ds[0] += 0.5 * std::abs(cell.nei[i]->nor.x);ds[1] += 0.5 * std::abs(cell.nei[i]->nor.y);}
    std::vector<double> LMDc = {(std::abs(cell.phy.u) + cell.phy.a) * ds[0],(std::abs(cell.phy.v) + cell.phy.a) * ds[1]};
    std::vector<double> LMDv = {(cc::gamma/cell.phy.rho) * (cell.phy.mu/cc::Pr + cell.tur.mut/SA::Prt) * (ds[0]*ds[0]/cell.vol),
                                (cc::gamma/cell.phy.rho) * (cell.phy.mu/cc::Pr + cell.tur.mut/SA::Prt) * (ds[1]*ds[1]/cell.vol)};
    cell.localdt =  fatime::CFL * cell.vol / (LMDc[0] + LMDc[1] + fatime::C * LMDv[0] + fatime::C * LMDv[1]); 
}