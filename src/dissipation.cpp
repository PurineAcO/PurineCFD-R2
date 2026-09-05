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
    for(int j=0;j<4;j++){cell.disspiation.L[j] = 0.0;}
    for(int i=0;i<cell.ecnt;i++){
        cc::cell_class* nei = cc::find_neicell(cell.index, cell.nei[i]);
        if(nei == NULL){continue;}
        for(int j=0;j<4;j++){
            cell.disspiation.L[j] += nei->conser[j] - cell.conser[j];
        }
    }
}

void jst::JST_dissipation(cc::cell_class &cell){
    for(int j=0;j<4;j++){cell.disspiation.Fd[j] = 0.0;}
    for(int i=0;i<cell.ecnt;i++){
        cc::face_class* face = cell.nei[i];
        cc::cell_class* nei = cc::find_neicell(cell.index, face);
        if(nei == NULL){continue;}

        // 该面谱半径(带面长): |u*nx+v*ny| + a*|n|
        double nx = face->nor.x, ny = face->nor.y;
        double nlen = std::hypot(nx, ny);
        double un = face->phy.u*nx + face->phy.v*ny;
        double lam = std::abs(un) + face->phy.a*nlen;
        if(lam < 1e-30){ continue; }

        double eps2 = jst::k2 * std::max(cell.disspiation.Y, nei->disspiation.Y);
        double eps4 = std::max(0.0, jst::k4 - eps2);

        for(int j=0;j<4;j++){
            cell.disspiation.Fd[j] += lam * eps2 * (nei->conser[j] - cell.conser[j]);
            cell.disspiation.Fd[j] += lam * eps4 * (cell.disspiation.L[j] - nei->disspiation.L[j]);
        }
    }
}
