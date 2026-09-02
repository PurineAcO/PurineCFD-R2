#include "dissipation.h"
#include "classconfig.h"
#include <cmath>

// JST: 二阶抑制激波, 四阶抑制锯齿
void jst::shockwave_recognize(cc::cell_class& cell){
    double up = 0, down = 0;
    for(int i=0;i<cell.ecnt;i++){
        cc::cell_class* nei = cc::find_neicell(cell.index, cell.nei[i]);
        if(nei == NULL){continue;}
        up += std::abs(nei->phy.p - cell.phy.p);
        down += nei->phy.p + cell.phy.p;
    }
    cell.disspiation.Y = (down > 0) ? up/down : 0.0;
}

void jst::laplace_dissipation(cc::cell_class &cell){
    for(int j=0;j<cc::NEQ;j++){cell.disspiation.L[j] = 0.0;}
    for(int i=0;i<cell.ecnt;i++){
        cc::cell_class* nei = cc::find_neicell(cell.index, cell.nei[i]);
        if(nei == NULL){continue;}
        for(int j=0;j<cc::NEQ;j++){
            cell.disspiation.L[j] += nei->conser[j] - cell.conser[j];
        }
    }
}

void jst::JST_dissipation(cc::cell_class &cell){
    for(int j=0;j<cc::NEQ;j++){cell.disspiation.Fd[j] = 0.0;}
    for(int i=0;i<cell.ecnt;i++){
        cc::cell_class* nei = cc::find_neicell(cell.index, cell.nei[i]);
        if(nei == NULL){continue;}

        double LMD_me = 0, LMD_nei = 0;
        for(int j=0;j<cell.ecnt;j++){
            LMD_me += std::abs(cell.nei[j]->phy.u * cell.nei[j]->nor.x + cell.nei[j]->phy.v * cell.nei[j]->nor.y);
            LMD_me += cell.nei[j]->phy.a * std::hypot(cell.nei[j]->nor.x, cell.nei[j]->nor.y);
        }
        for(int j=0;j<nei->ecnt;j++){
            LMD_nei += std::abs(nei->nei[j]->phy.u * nei->nei[j]->nor.x + nei->nei[j]->phy.v * nei->nei[j]->nor.y);
            LMD_nei += nei->nei[j]->phy.a * std::hypot(nei->nei[j]->nor.x, nei->nei[j]->nor.y);
        }
        double LMD = 0.5 * (LMD_me + LMD_nei);

        double eps2 = jst::k2 * std::max(cell.disspiation.Y, nei->disspiation.Y);
        double eps4 = std::max(0.0, jst::k4 - eps2);

        for(int j=0;j<cc::NEQ;j++){
            cell.disspiation.Fd[j] += LMD * eps2 * (nei->conser[j] - cell.conser[j]);
            cell.disspiation.Fd[j] += LMD * eps4 * (cell.disspiation.L[j] - nei->disspiation.L[j]);
        }
    }
}
