#pragma once

#include "classconfig.hpp"
#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

namespace res {
    // 流动方程各量相对首次检查的下降比例低于该值即认为收敛(不含湍流工作变量)
    inline constexpr double drop_target = 1e-4;
    // 报告各物理量的全场最大变化量及位置
    void report_update(int step);
    // rho,u,v,e 相对首次检查下降比例中最差的一个, 越小越收敛
    double worst_drop();
}

namespace res {

static constexpr int fields = 5; // rho,u,v,e,miubl
static std::vector<std::array<double,fields>> former; // 上次检查时的状态
static bool first = true;
static double max_update[fields] = {};
static double first_update[fields] = {}; // 首次检查时的最大变化量
static bool have_first = false;

inline void report_update(int step){
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

inline double worst_drop(){
    // 收敛判据只看流动方程(rho,u,v,e), 湍流工作变量 miubl 不参与
    constexpr int checked = fields - 1;
    if(!have_first){
        for(int s=0;s<checked;s++){
            first_update[s] = max_update[s];
        }
        have_first = true;
        return 1.0;
    }
    double worst = 0.0;
    for(int s=0;s<checked;s++){
        worst = std::max(worst,max_update[s]/std::max(first_update[s],1e-300));
    }
    return worst;
}

}
