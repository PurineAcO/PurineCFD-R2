#pragma once

#include "classconfig.hpp"
#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

namespace res {
    // 流动方程各量相对首次检查的下降比例低于该值即认为收敛(不含湍流工作变量)
    inline constexpr double drop_target = 1e-6;
    // 报告各物理量的L2变化量(全场变化量平方和除以单元数再开方)及变化最大单元
    void report_update(int step);
    // rho,u,v,e 相对首次检查下降比例中最差的一个, 越小越收敛
    double worst_drop();
}

namespace res {

static constexpr int fields = 5; // rho,u,v,e,miubl
static std::vector<std::array<double,fields>> former; // 上次检查时的状态
static bool first = true;
static double l2_update[fields] = {};
static double first_update[fields] = {}; // 首次检查时的L2变化量
static bool have_first = false;

inline void report_update(int step){
    if(former.empty()){
        former.resize(cc::cell_num);
    }
    double sum_square[fields] = {};
    for(int i=0;i<cc::cell_num;i++){
        cc::cell_class& cell = cc::CellList[i];
        double state[fields] = {cell.phy.rho,cell.phy.u,cell.phy.v,cell.otphy.e,cell.tur.miubl};
        if(first){
            for(int s=0;s<fields;s++){
                former[i][s] = state[s];
            }
            continue;
        }
        for(int s=0;s<fields;s++){
            const double change = state[s] - former[i][s];
            sum_square[s] += change*change;
            former[i][s] = state[s];
        }
    }
    if(first){
        first = false;
        return;
    }
    for(int s=0;s<fields;s++){
        l2_update[s] = std::sqrt(sum_square[s]/cc::cell_num);
    }
    static bool header = false;
    if(!header){
        printf("%7s %12s %12s %12s %12s %12s\n","step","drho","du","dv","de","dnu_tilde");
        header = true;
    }
    printf("%7d",step);
    for(int s=0;s<fields;s++){
        printf(" %12.6e",l2_update[s]);
    }
    printf("\n");
}

inline double worst_drop(){
    // 收敛判据只看流动方程(rho,u,v,e), 湍流工作变量 miubl 不参与
    constexpr int checked = fields - 1;
    if(!have_first){
        for(int s=0;s<checked;s++){
            first_update[s] = l2_update[s];
        }
        have_first = true;
        return 1.0;
    }
    double worst = 0.0;
    for(int s=0;s<checked;s++){
        worst = std::max(worst,l2_update[s]/std::max(first_update[s],1e-300));
    }
    return worst;
}

}
