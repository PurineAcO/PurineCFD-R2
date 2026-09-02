#include "residual.h"
#include "classconfig.h"
#include "config.h"
#include <cmath>
#include <array>
#include <vector>
#include <cstdio>

namespace res {

namespace {
    std::vector<std::array<double,cc::NEQ>> g_prev;
    bool g_first = true;
}

void report_update(int step){
    if(g_prev.empty()){ g_prev.resize(cc::cell_num); }
    double maxd[cc::NEQ] = {0,0,0,0};
    double maxnorm = -1.0;
    int maxcell = 0;
    double maxx = 0.0, maxy = 0.0;
    std::size_t idx = 0;
    for(cc::cell_class& cell : cc::CellList){
        double prim[cc::NEQ] = {cell.phy.rho, cell.phy.u, cell.phy.v, cell.phy.e};
        if(g_first){
            for(int s=0;s<cc::NEQ;s++){ g_prev[idx][s] = prim[s]; }
        }else{
            double d[cc::NEQ]; double rn = 0.0;
            for(int s=0;s<cc::NEQ;s++){ d[s] = prim[s] - g_prev[idx][s]; rn += d[s]*d[s]; }
            rn = std::sqrt(rn);
            if(rn > maxnorm){
                maxnorm = rn;
                for(int s=0;s<cc::NEQ;s++){ maxd[s] = d[s]; }
                maxcell = cell.index; maxx = cell.center.x; maxy = cell.center.y;
            }
            for(int s=0;s<cc::NEQ;s++){ g_prev[idx][s] = prim[s]; }
        }
        idx++;
    }
    if(g_first){ g_first = false; return; }

    static bool header = false;
    if(!header){
        printf("%6s %12s %12s %12s %12s  %s\n",
               "step", "drho", "du", "dv", "de", "max残差位置");
        header = true;
    }
    printf("%6d", step);
    for(int s=0;s<cc::NEQ;s++){
        printf(" %12.6e", maxd[s]);
    }
    printf("  cell#%d (x=%.3f,y=%.6f)\n", maxcell, maxx, maxy);
}

}
