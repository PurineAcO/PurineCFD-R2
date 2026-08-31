#include "classconfig.h"
#include "config.h"


void velocity_inlet_boundary(){
    for(cc::face_class* face : cc::VILFaces){
        face->phy.u = cc::VIL_DEFINE.u;
        face->phy.v = cc::VIL_DEFINE.v;
        face->phy.T = cc::VIL_DEFINE.T;
        face->tur.miubl = cc::VIL_DEFINE.tur.miubl;
        face->phy.p = face->nei[0] ? face->nei[0]->phy.p : face->nei[1]->phy.p;
        face->phy.rho = face->phy.p / cc::R / face->phy.T;  // 强行满足本构关系
        face->phy.rhograd = face->tur.miublgrad = {0.0,0.0};
    }
}

void wall_boundary(){
    for(cc::face_class* wall : cc::WallFaces){
        wall->phy.u = 0;wall->phy.v = 0;wall->tur.miubl = 0;
        wall->phy.T = (wall->nei[0] != NULL) ? wall->nei[0]->phy.T : wall->nei[1]->phy.T;
        wall->phy.rho = (wall->nei[0] != NULL) ? wall->nei[0]->phy.rho : wall->nei[1]->phy.rho;
        wall->phy.Tgrad = wall->phy.rhograd = {0.0,0.0};
    }
}

void pressure_outlet_boundary(){
    for(cc::face_class *pol : cc::POLFaces){
        cc::cell_class *cell_in = pol->nei[0] ? pol->nei[0] : pol->nei[1];
        pol->phy.u = cell_in->phy.u;
        pol->phy.v = cell_in->phy.v;
        pol->phy.p = cc::POL_DEFINE.p;
        pol->phy.ugrad = pol->phy.vgrad = pol->phy.Tgrad = pol->phy.rhograd = pol->tur.miublgrad = {0.0,0.0};
        bool outer;for(int i=0;i<cell_in->ecnt;i++){if(cell_in->nei[i] == pol){outer = cell_in->fnorm[i];}}
        if((2*outer-1)*cc::dot(cc::vec2{pol->phy.u, pol->phy.v}, pol->nor) > 0){// 正常流出
            pol->phy.rho = cell_in->phy.rho;
            pol->phy.T = pol->phy.p/cc::R/pol->phy.rho;
            pol->tur.miubl = cell_in->tur.miubl;}
        else{// 回流
            pol->phy.T = cc::POL_DEFINE.T - (pol->phy.u*pol->phy.u + pol->phy.v*pol->phy.v) / (2*cc::Cp);
            pol->phy.rho = pol->phy.p / cc::R / pol->phy.T;
            pol->tur.miubl = cc::POL_DEFINE.tur.miubl;}
    }
}