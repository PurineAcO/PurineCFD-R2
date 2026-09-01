#include "dissipation.h"
#include "classconfig.h"
#include <vector>
#include <cmath>

inline double theta_(cc::cell_class* me,cc::cell_class* nei){
    double t = me->disspiation.lambda.x * (nei->center.x - me->center.x)
             + me->disspiation.lambda.y * (nei->center.y - me->center.y);
    if(t < 0.0){ t = 0.0; }
    if(t > 2.0){ t = 2.0; }
    return t;
}

void jst::JST_dissipation_INIT(cc::cell_class &cell){
    for(int i=0;i<cell.ecnt;i++){
        cc::cell_class* nei = cc::find_neicell(cell.index, cell.nei[i]);
        if(nei == NULL){continue;}
        cell.disspiation.I[0] += (nei->center.x - cell.center.x)*(nei->center.x - cell.center.x);
        cell.disspiation.I[1] += (nei->center.x - cell.center.x)*(nei->center.y - cell.center.y);
        cell.disspiation.I[2] += (nei->center.y - cell.center.y)*(nei->center.y - cell.center.y);
        cell.disspiation.R[0] += nei->center.x - cell.center.x;
        cell.disspiation.R[1] += nei->center.y - cell.center.y;
    }
    double d = cell.disspiation.I[0]*cell.disspiation.I[2] - cell.disspiation.I[1]*cell.disspiation.I[1];
    cell.disspiation.lambda.x = (-1) * (cell.disspiation.R[0]*cell.disspiation.I[2]-cell.disspiation.R[1]*cell.disspiation.I[1]) / d;
    cell.disspiation.lambda.y = (-1) * (cell.disspiation.R[1]*cell.disspiation.I[0]-cell.disspiation.R[0]*cell.disspiation.I[1]) / d;
}

void jst::shockwave_recognize(cc::cell_class& cell){
    double up = 0;double down = 0;
    for(int i=0;i<cell.ecnt;i++){
        cc::cell_class* nei = cc::find_neicell(cell.index, cell.nei[i]);
        if(nei == NULL){continue;}
        double theta = theta_(&cell, nei);
        up += theta * (nei->phy.p - cell.phy.p);
        down += nei->phy.p + cell.phy.p;
    }
    cell.disspiation.Y = up/down;
}

void jst::laplace_dissipation(cc::cell_class &cell){
    cell.disspiation.L = {0.0,0.0,0.0,0.0,0.0};
    for(int i=0;i<cell.ecnt;i++){
        cc::cell_class* nei = cc::find_neicell(cell.index, cell.nei[i]);
        if(nei == NULL){continue;}
        double theta = theta_(&cell,nei);
        for(int j=0;j<5;j++){
            cell.disspiation.L[j] += theta * (nei->conser[j] - cell.conser[j]);
        }
    }
}

void jst::JST_dissipation(cc::cell_class &cell){
    cell.disspiation.Fd = {0.0,0.0,0.0,0.0,0.0};
    for(int i=0;i<cell.ecnt;i++){
        cc::cell_class* nei = cc::find_neicell(cell.index, cell.nei[i]);
        if(nei == NULL){continue;}

        // 计算几何权
        double theta = theta_(&cell, nei);

        // 建立谱半径 (对流通量特征值上界 |u·n| + a|n|)
        double LMD_me = 0;double LMD_nei = 0;double LMD = 0;
        for(int j=0;j<cell.ecnt;j++){ // 本家网格
        LMD_me += std::abs(cell.nei[j]->phy.u * cell.nei[j]->nor.x + cell.nei[j]->phy.v * cell.nei[j]->nor.y);
        LMD_me += cell.nei[j]->phy.a * sqrt(cell.nei[j]->nor.x*cell.nei[j]->nor.x + cell.nei[j]->nor.y*cell.nei[j]->nor.y);}
        for(int j=0;j<nei->ecnt;j++){ // 邻接网格
        LMD_nei += std::abs(nei->nei[j]->phy.u * nei->nei[j]->nor.x + nei->nei[j]->phy.v * nei->nei[j]->nor.y);
        LMD_nei += nei->nei[j]->phy.a * sqrt(nei->nei[j]->nor.x*nei->nei[j]->nor.x + nei->nei[j]->nor.y*nei->nei[j]->nor.y);}
        LMD = 0.5 * (LMD_me + LMD_nei);

        // 建立耗散系数
        double eps2 = jst::k2 * std::max(cell.disspiation.Y,nei->disspiation.Y);
        double eps4 = std::max(0.0,(jst::k4-eps2));

        // 形成二阶项和四阶项
        for(int j=0;j<5;j++){
            cell.disspiation.Fd[j] += LMD * eps2 * theta * (nei->conser[j] - cell.conser[j]);
            cell.disspiation.Fd[j] += LMD * eps4 * theta * (nei->disspiation.L[j] - cell.disspiation.L[j]);
        }
    }
}