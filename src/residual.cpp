#include "residual.h"
#include "classconfig.h"
#include "physic.h"
#include "udf.h"
#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

namespace res {

namespace {
    constexpr int fields = 5; // rho,u,v,e,miubl
    std::vector<std::array<double,fields>> former; // 上次检查时的状态
    bool first = true;
    double max_update[fields] = {};
}

void report_update(int step){
    if(former.empty()){
        former.resize(cc::cell_num);
    }
    double delta[fields] = {0,0,0,0,0};
    double worst_norm = -1.0;
    int worst_cell = 0;
    for(int i=0;i<cc::cell_num;i++){
        cc::cell_class& cell = cc::CellList[i];
        double state[fields] = {cell.phy.rho,cell.phy.u,cell.phy.v,cell.phy.e,cell.tur.miubl};
        if(first){
            for(int s=0;s<fields;s++){
                former[i][s] = state[s];
            }
            continue;
        }
        double norm = 0.0;
        for(int s=0;s<fields;s++){
            double change = state[s] - former[i][s];
            norm += change*change;
            if(std::fabs(change) > delta[s]){
                delta[s] = std::fabs(change);
            }
            former[i][s] = state[s];
        }
        norm = std::sqrt(norm);
        if(norm > worst_norm){
            worst_norm = norm;
            worst_cell = cell.index;
        }
    }
    if(first){
        first = false;
        return;
    }
    for(int s=0;s<fields;s++){
        max_update[s] = delta[s];
    }
    static bool header = false;
    if(!header){
        printf("%7s %12s %12s %12s %12s %12s  %s\n","step","drho","du","dv","de","dnu_tilde",
               "maxcell");
        header = true;
    }
    printf("%7d",step);
    for(int s=0;s<fields;s++){
        printf(" %12.6e",max_update[s]);
    }
    cc::cell_class& cell = cc::gotocell(worst_cell);
    printf("  #%d(%.4f,%.4f)\n",worst_cell,cell.center.x,cell.center.y);
}

double relative_update(){
    double rho = FAR_DEFINE.p/(cc::R*FAR_DEFINE.T);
    double U = std::max(std::hypot(FAR_DEFINE.u,FAR_DEFINE.v),1e-20);
    double scale[fields] = {rho,U,U,cc::Cv*FAR_DEFINE.T,
                            3.0*sutherland::dynamic_viscosity(FAR_DEFINE.T)/rho};
    double maximum = 0.0;
    for(int s=0;s<fields;s++){
        maximum = std::max(maximum,max_update[s]/scale[s]);
    }
    return maximum;
}

}
