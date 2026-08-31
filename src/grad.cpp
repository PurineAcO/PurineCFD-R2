#include "classconfig.h"

void green_gauss_cell_based(cc::cell_class & cell){
    cell.phy.rhograd    = {0.0, 0.0};
    cell.phy.ugrad      = {0.0, 0.0};
    cell.phy.vgrad      = {0.0, 0.0};
    cell.phy.Tgrad      = {0.0, 0.0};
    cell.tur.miublgrad  = {0.0, 0.0};
    for(int i=0;i<cell.ecnt;i++){
        // cell.nei[i]->face_physic_mid();  如果采用高阶重构,此处不该是有东西的.
        const double s = (2*cell.fnorm[i]-1) / cell.vol;
        cell.phy.rhograd   += (cell.nei[i]->phy.rho   * s) * cell.nei[i]->nor;
        cell.phy.ugrad     += (cell.nei[i]->phy.u     * s) * cell.nei[i]->nor;
        cell.phy.vgrad     += (cell.nei[i]->phy.v     * s) * cell.nei[i]->nor;
        cell.phy.Tgrad     += (cell.nei[i]->phy.T     * s) * cell.nei[i]->nor;
        cell.tur.miublgrad += (cell.nei[i]->tur.miubl * s) * cell.nei[i]->nor;
    }
}