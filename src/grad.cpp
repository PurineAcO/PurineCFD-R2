#include "classconfig.h"

void green_gauss_cell_based_physics(cc::cell_class & cell){
    cell.phy.rhograd    = {0.0, 0.0};
    cell.phy.ugrad      = {0.0, 0.0};
    cell.phy.vgrad      = {0.0, 0.0};
    cell.phy.Tgrad      = {0.0, 0.0};
    cell.tur.miublgrad  = {0.0, 0.0};
    for(int i=0;i<cell.ecnt;i++){
        // cell.nei[i]->face_physic_mid();  如果采用高阶重构,此处不该是有东西的.
        cell.phy.rhograd.x += cell.nei[i]->phy.rho * (2*cell.fnorm[i]-1) *cell.nei[i]->nor.x / cell.vol;
        cell.phy.rhograd.y += cell.nei[i]->phy.rho * (2*cell.fnorm[i]-1) *cell.nei[i]->nor.y / cell.vol;
        cell.phy.ugrad.x += cell.nei[i]->phy.u * (2*cell.fnorm[i]-1) *cell.nei[i]->nor.x / cell.vol;
        cell.phy.ugrad.y += cell.nei[i]->phy.u * (2*cell.fnorm[i]-1) *cell.nei[i]->nor.y / cell.vol;
        cell.phy.vgrad.x += cell.nei[i]->phy.v * (2*cell.fnorm[i]-1) *cell.nei[i]->nor.x / cell.vol;
        cell.phy.vgrad.y += cell.nei[i]->phy.v * (2*cell.fnorm[i]-1) *cell.nei[i]->nor.y / cell.vol;
        cell.phy.Tgrad.x += cell.nei[i]->phy.T * (2*cell.fnorm[i]-1) *cell.nei[i]->nor.x / cell.vol;
        cell.phy.Tgrad.y += cell.nei[i]->phy.T * (2*cell.fnorm[i]-1) *cell.nei[i]->nor.y / cell.vol;
        cell.tur.miublgrad.x += cell.nei[i]->tur.miubl * (2*cell.fnorm[i]-1) *cell.nei[i]->nor.x / cell.vol;
        cell.tur.miublgrad.y += cell.nei[i]->tur.miubl * (2*cell.fnorm[i]-1) *cell.nei[i]->nor.y / cell.vol;
    }
}