#include "residual.h"
#include "classconfig.h"
#include "config.h"
#include <cmath>
#include <array>
#include <vector>
#include <cstdio>

namespace res {

namespace {
    std::vector<std::array<double,4>> g_prev;
    bool g_first = true;
    double g_rms = 1e30;
}

// 全场场更新RMS(密度), 用于定常收敛判据
double current_residual(){ return g_rms; }

void report_update(int step){
    if(g_prev.empty()){ g_prev.resize(cc::cell_num); }
    double maxd[4] = {0,0,0,0};
    double maxnorm = -1.0;
    int maxcell = 0;
    double maxx = 0.0, maxy = 0.0;
    double sum = 0.0; int cnt = 0;
    std::size_t idx = 0;
    for(cc::cell_class& cell : cc::CellList){
        double prim[4] = {cell.phy.rho, cell.phy.u, cell.phy.v, cell.phy.e};
        if(g_first){
            for(int s=0;s<4;s++){ g_prev[idx][s] = prim[s]; }
        }else{
            double d[4]; double rn = 0.0;
            for(int s=0;s<4;s++){ d[s] = prim[s] - g_prev[idx][s]; rn += d[s]*d[s]; }
            rn = std::sqrt(rn);
            sum += d[0]*d[0]; cnt++;
            if(rn > maxnorm){
                maxnorm = rn;
                for(int s=0;s<4;s++){ maxd[s] = d[s]; }
                maxcell = cell.index; maxx = cell.center.x; maxy = cell.center.y;
            }
            for(int s=0;s<4;s++){ g_prev[idx][s] = prim[s]; }
        }
        idx++;
    }
    if(g_first){ g_first = false; return; }
    g_rms = (cnt > 0) ? std::sqrt(sum/cnt) : 1e30;

    static bool header = false;
    if(!header){
        printf("%6s %12s %12s %12s %12s  %s\n",
               "step", "drho", "du", "dv", "de", "max残差位置");
        header = true;
    }
    printf("%6d", step);
    for(int s=0;s<4;s++){
        printf(" %12.6e", maxd[s]);
    }
    printf("  cell#%d (x=%.3f,y=%.6f)\n", maxcell, maxx, maxy);
}

}
