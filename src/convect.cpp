#include "classconfig.h"
#include "config.h"
#include <vector>

void convect_JST(cc::cell_class &cell){
    cell.convect = {0.0,0.0,0.0,0.0,0.0};
    for(int i=0;i<cell.ecnt;i++){
        cc::face_class* face = cell.nei[i];short outer = 2*cell.fnorm[i]-1;
        std::vector<double> F = {face->phy.rho*face->phy.u,
                                face->phy.rho*face->phy.u*face->phy.u + face->phy.rho*cc::R*face->phy.T,
                                face->phy.rho*face->phy.u*face->phy.v,
                                face->phy.rho*face->phy.u*(cc::Cp*face->phy.T + 0.5*(face->phy.u*face->phy.u+face->phy.v*face->phy.v)),
                                face->phy.rho*face->phy.u*face->tur.miubl};
        std::vector<double> G = {
                                face->phy.rho*face->phy.v,
                                face->phy.rho*face->phy.u*face->phy.v,
                                face->phy.rho*face->phy.v*face->phy.v + face->phy.rho*cc::R*face->phy.T,
                                face->phy.rho*face->phy.v*(cc::Cp*face->phy.T + 0.5*(face->phy.u*face->phy.u+face->phy.v*face->phy.v)),
                                face->phy.rho*face->phy.v*face->tur.miubl};
        for(int j=0;j<=4;j++){cell.convect[j] += F[j]*outer*face->nor.x + G[j]*outer*face->nor.y;}
    }
}