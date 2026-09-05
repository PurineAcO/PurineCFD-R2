#include "classconfig.h"
#include "config.h"

// 无粘通量
void convect_JST(cc::cell_class &cell){
    for(int j=0;j<4;j++){cell.convect[j] = 0.0;}
    for(int i=0;i<cell.ecnt;i++){
        cc::face_class* face = cell.nei[i];
        short outer = 2*cell.fnorm[i]-1;
        double rho = face->phy.rho, u = face->phy.u, v = face->phy.v;
        double H = cc::Cp*face->phy.T + 0.5*(u*u+v*v); 
        double p = cc::R*face->phy.rho*face->phy.T;
        double F[4] = {rho*u, rho*u*u + p, rho*u*v, rho*u*H};
        double G[4] = {rho*v, rho*u*v, rho*v*v + p, rho*v*H};
        for(int j=0;j<4;j++){cell.convect[j] += outer * (F[j]*face->nor.x + G[j]*face->nor.y);}
    }
}
