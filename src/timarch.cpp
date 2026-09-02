#include "classconfig.h"
#include "config.h"
#include "timarch.h"
#include <cmath>

void local_timestep(cc::cell_class &cell){
    double dsx = 0, dsy = 0;
    for(int i=0;i<cell.ecnt;i++){
        dsx += 0.5 * std::abs(cell.nei[i]->nor.x);
        dsy += 0.5 * std::abs(cell.nei[i]->nor.y);
    }
    double LMDx = (std::abs(cell.phy.u) + cell.phy.a) * dsx;
    double LMDy = (std::abs(cell.phy.v) + cell.phy.a) * dsy;
    cell.localdt = fatime::CFL * cell.vol / (LMDx + LMDy);
}
