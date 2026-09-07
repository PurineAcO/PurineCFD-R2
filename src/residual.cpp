#include "residual.h"
#include "classconfig.h"
#include "config.h"
#include <cmath>
#include <array>
#include <vector>
#include <cstdio>

namespace res {

namespace {
    constexpr int NP = 5;   // rho,u,v,e,miubl
    std::vector<std::array<double,NP>> g_prev;
    bool g_first = true;
    double g_m2[NP];    // 均值平方残差
    double g_amax[NP];  // 每物理量最大绝对值残差
}

double current_residual_all(){
    double m = g_m2[0];
    for(int s=1;s<4;s++){ if(g_m2[s] > m) m = g_m2[s]; }
    return m;
}

double residual_absmax(){
    double m = g_amax[0];
    for(int s=1;s<4;s++){ if(g_amax[s] > m) m = g_amax[s]; }
    return m;
}

void report_update(int step){
    if(g_prev.empty()){ g_prev.resize(cc::cell_num); }
    double maxd[NP] = {0,0,0,0,0};
    double maxnorm = -1.0;
    int maxcell = 0;
    double maxx = 0.0, maxy = 0.0;
    double sumSq[NP] = {0,0,0,0,0};
    int cnt = 0;
    std::size_t idx = 0;
    for(cc::cell_class& cell : cc::CellList){
        double prim[NP] = {cell.phy.rho, cell.phy.u, cell.phy.v, cell.phy.e, cell.tur.miubl};
        if(g_first){
            for(int s=0;s<NP;s++){ g_prev[idx][s] = prim[s]; }
        }else{
            double d[NP]; double rn = 0.0;
            for(int s=0;s<NP;s++){ d[s] = prim[s] - g_prev[idx][s]; rn += d[s]*d[s]; }
            rn = std::sqrt(rn);
            for(int s=0;s<NP;s++){ sumSq[s] += d[s]*d[s]; }
            for(int s=0;s<NP;s++){ if(std::fabs(d[s]) > maxd[s]){ maxd[s] = std::fabs(d[s]); } }
            cnt++;
            if(rn > maxnorm){
                maxnorm = rn;
                maxcell = cell.index; maxx = cell.center.x; maxy = cell.center.y;
            }
            for(int s=0;s<NP;s++){ g_prev[idx][s] = prim[s]; }
        }
        idx++;
    }
    if(g_first){ g_first = false; return; }
    if(cnt > 0){
        for(int s=0;s<NP;s++){ g_m2[s] = sumSq[s]/cnt; g_amax[s] = maxd[s]; }
    }

    static bool header = false;
    if(!header){
        printf("%7s %12s %12s %12s %12s %12s  %s\n",
               "step", "drho", "du", "dv", "de", "dmiubl", "maxcell");
        header = true;
    }
    printf("%7d", step);
    for(int s=0;s<NP;s++){ printf(" %12.6e", g_amax[s]); }
    printf("  #%d(%.4f,%.4f)\n", maxcell, maxx, maxy);
}

}
