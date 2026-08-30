#include "classconfig.h"
#include "config.h"


void velocity_inlet_boundary(){
    for(cc::face_class* face : cc::VILFaces){
        face->phy.u = cc::VIL_DEFINE.u;
        face->phy.v = cc::VIL_DEFINE.v;
        face->phy.T = cc::VIL_DEFINE.T;
        face->tur.miubl = cc::VIL_DEFINE.tur.miubl;
        face->phy.p = face->nei[0] ? face->nei[0]->phy.p : face->nei[1]->phy.p;
    }
}