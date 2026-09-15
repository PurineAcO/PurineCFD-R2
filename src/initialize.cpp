#include "initialize.h"
#include "classconfig.h"
#include "physic.h"
#include "udf.h"
#include <cstdio>

void std_initialize(){
    double rho_inf = FAR_DEFINE.p/(cc::R*FAR_DEFINE.T);
    double miubl_inf = 3.0*sutherland::dynamic_viscosity(FAR_DEFINE.T)/rho_inf;
    for(cc::cell_class& cell : cc::CellList){
        cell.phy.T = FAR_DEFINE.T;
        cell.phy.p = FAR_DEFINE.p;
        cell.phy.u = FAR_DEFINE.u;
        cell.phy.v = FAR_DEFINE.v;
        cell.phy.rho = cell.phy.p/cell.phy.T/cc::R;
        cell.phy.a = get_sonic_velocity(cell.phy.T);
        cell.phy.e = cc::Cv*cell.phy.T +
                     0.5*(cell.phy.u*cell.phy.u + cell.phy.v*cell.phy.v);
        cell.tur.miubl = miubl_inf;
    }
    for(cc::face_class& face : cc::FaceList){
        face.tur.miubl = (face.type == cc::WALL) ? 0.0 : miubl_inf;
        if(face.type == cc::INTER){
            face.face_physic_mid();
        }
    }
    printf("STD Initialization OK!, u is: %f\n",cc::CellList[0].phy.u);
}
