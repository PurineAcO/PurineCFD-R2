#include "residual.h"
#include "classconfig.h"
#include "config.h"
#include <cmath>
#include <array>
#include <vector>
#include <cstdio>

namespace res {

namespace {
    std::vector<std::array<double,5>> g_prev;
    bool g_first = true;
}

void report_update(int step){
    if(g_prev.empty()){ g_prev.resize(cc::cell_num); }
    double resSq[5] = {0,0,0,0,0};
    std::size_t idx = 0;
    for(cc::cell_class& cell : cc::CellList){
        double prim[5] = {cell.phy.rho, cell.phy.u, cell.phy.v, cell.phy.e, cell.tur.miubl};
        if(g_first){
            for(int s=0;s<5;s++){ g_prev[idx][s] = prim[s]; }
        }else{
            for(int s=0;s<5;s++){
                double d = prim[s] - g_prev[idx][s];
                resSq[s] += d*d;
            }
            for(int s=0;s<5;s++){ g_prev[idx][s] = prim[s]; }
        }
        idx++;
    }
    if(g_first){ g_first = false; return; }

    static bool header = false;
    if(!header){
        printf("%6s %12s %12s %12s %12s %12s\n",
               "step", "drho", "du", "dv", "de", "dnu");
        header = true;
    }
    double N = (double)cc::cell_num;
    printf("%6d", step);
    for(int s=0;s<5;s++){
        printf(" %12.6e", std::sqrt(resSq[s]/N));
    }
    printf("\n");
}

}
