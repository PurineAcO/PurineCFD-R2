#include "convect.h"

void convect_JST(cc::face_class& face){
    const double rho = face.phy.rho,u = face.phy.u,v = face.phy.v;
    const double H = cc::Cp*face.phy.T + 0.5*(u*u + v*v);
    const double p = cc::R*rho*face.phy.T;
    const double F[4] = {rho*u,rho*u*u + p,rho*u*v,rho*u*H};
    const double G[4] = {rho*v,rho*u*v,rho*v*v + p,rho*v*H};
    for(int j=0;j<4;j++){
        face.convect[j] = F[j]*face.nor.x + G[j]*face.nor.y;
    }
}

void assemble_flux(cc::cell_class& cell){
    for(int j=0;j<4;j++){
        cell.convect[j] = 0.0;
        cell.visflux[j] = 0.0;
    }
    for(int i=0;i<cell.ecnt;i++){
        cc::face_class& face = *cell.faces[i];
        const int outer = 2*cell.fnorm[i] - 1;
        for(int j=0;j<4;j++){
            cell.convect[j] += outer*face.convect[j];
            cell.visflux[j] += outer*face.visflux[j];
        }
    }
}
