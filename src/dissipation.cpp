#include "dissipation.h"
#include <cmath>

void jst::shockwave_recognize(cc::cell_class& cell){
    double up = 0,down = 0;
    for(int i=0;i<cell.ecnt;i++){
        cc::cell_class* neighbor = cell.nei[i];
        if(neighbor == nullptr){
            continue;
        }
        up += std::abs(neighbor->phy.p - cell.phy.p);
        down += neighbor->phy.p + cell.phy.p;
    }
    cell.diss.Y = (down > 0) ? up/down : 0.0;
}

void jst::laplace_dissipation(cc::cell_class& cell){
    for(int j=0;j<4;j++){
        cell.diss.L[j] = 0.0;
    }
    for(int i=0;i<cell.ecnt;i++){
        cc::cell_class* neighbor = cell.nei[i];
        if(neighbor == nullptr){
            continue;
        }
        for(int j=0;j<4;j++){
            cell.diss.L[j] += neighbor->conser[j] - cell.conser[j];
        }
    }
}

void jst::JST_dissipation(cc::cell_class& cell){
    for(int j=0;j<4;j++){
        cell.diss.Fd[j] = 0.0;
    }
    for(int i=0;i<cell.ecnt;i++){
        cc::face_class* face = cell.faces[i];
        cc::cell_class* neighbor = cell.nei[i];
        if(neighbor == nullptr){
            continue;
        }
        const double lam = face->lam;
        if(lam < 1e-30){
            continue;
        }
        double eps2 = jst::k2*std::max(cell.diss.Y,neighbor->diss.Y);
        double eps4 = std::max(0.0,jst::k4 - eps2);
        for(int j=0;j<4;j++){
            cell.diss.Fd[j] += lam*eps2*(neighbor->conser[j] - cell.conser[j]);
            cell.diss.Fd[j] += lam*eps4*(cell.diss.L[j] - neighbor->diss.L[j]);
        }
    }
}
