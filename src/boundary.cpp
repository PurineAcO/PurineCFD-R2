#include "classconfig.h"
#include "config.h"

void velocity_inlet_boundary(){
    for(cc::face_class* face : cc::VILFaces){
        cc::cell_class* c = face->nei[0] ? face->nei[0] : face->nei[1];
        face->phy.u = cc::VIL_DEFINE.u;
        face->phy.v = cc::VIL_DEFINE.v;
        face->phy.T = cc::VIL_DEFINE.T;
        face->phy.p = c->phy.p;
        face->phy.rho = face->phy.p / cc::R / face->phy.T;
    }
}

void slip_wall_boundary(){
    for(cc::face_class* wall : cc::WallFaces){
        cc::cell_class* c = wall->nei[0] ? wall->nei[0] : wall->nei[1];
        double nx = wall->nor.x, ny = wall->nor.y;
        if(nx*(wall->mid.x - c->center.x) + ny*(wall->mid.y - c->center.y) < 0){ nx = -nx; ny = -ny; }
        double vn = (c->phy.u*nx + c->phy.v*ny) / (nx*nx + ny*ny);
        wall->phy.u = c->phy.u - vn*nx;
        wall->phy.v = c->phy.v - vn*ny;
        wall->phy.T = c->phy.T;
        wall->phy.rho = c->phy.rho;
    }
}

void pressure_outlet_boundary(){
    for(cc::face_class* pol : cc::POLFaces){
        cc::cell_class* c = pol->nei[0] ? pol->nei[0] : pol->nei[1];
        pol->phy.u = c->phy.u;
        pol->phy.v = c->phy.v;
        pol->phy.p = cc::POL_DEFINE.p;
        bool outer;for(int i=0;i<c->ecnt;i++){if(c->nei[i] == pol){outer = c->fnorm[i];}}
        if((2*outer-1)*cc::dot(cc::vec2{pol->phy.u, pol->phy.v}, pol->nor) > 0){
            pol->phy.rho = c->phy.rho;
            pol->phy.T = pol->phy.p/cc::R/pol->phy.rho;
        }else{
            pol->phy.T = cc::POL_DEFINE.T - (pol->phy.u*pol->phy.u + pol->phy.v*pol->phy.v) / (2*cc::Cp);
            pol->phy.rho = pol->phy.p / cc::R / pol->phy.T;
        }
    }
}
