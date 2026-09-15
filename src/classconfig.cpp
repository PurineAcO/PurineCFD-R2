#include "classconfig.h"
#include "physic.h"
#include <cmath>

namespace cc {

node_class::node_class(int number_,double x_,double y_):number(number_),x(x_),y(y_){}

face_class::face_class(int index_,int p1_,int p2_,int c1_,int c2_,short type_)
    :index(index_),type(type_),cell_1(c1_),cell_2(c2_){
    node[0] = &NodeList[p1_-1];
    node[1] = &NodeList[p2_-1];
    mid = vec2{0.5*(node[0]->x + node[1]->x),0.5*(node[0]->y + node[1]->y)};
    nor = vec2{node[0]->y - node[1]->y,node[1]->x - node[0]->x};
    len = std::hypot(nor.x,nor.y);
}

cell_class::cell_class(int index_,int f1_,int f2_,int f3_,int f4_)
    :index(index_),face{f1_,f2_,f3_,f4_}{}

void face_class::face_physic_mid(){
    if(type != INTER){
        return;
    }
    phy.rho = 0.5*(nei[0]->phy.rho + nei[1]->phy.rho);
    phy.u   = 0.5*(nei[0]->phy.u   + nei[1]->phy.u);
    phy.v   = 0.5*(nei[0]->phy.v   + nei[1]->phy.v);
    phy.T   = 0.5*(nei[0]->phy.T   + nei[1]->phy.T);
    tur.miubl = 0.5*(nei[0]->tur.miubl + nei[1]->tur.miubl);
}

void face_class::form_physic(){
    phy.e = get_energy(phy);
    phy.p = cc::R*phy.rho*phy.T;
    phy.a = get_sonic_velocity(phy.T);
}

void cell_class::form_physic(){
    phy.e = get_energy(phy);
    phy.p = cc::R*phy.rho*phy.T;
    phy.a = get_sonic_velocity(phy.T);
}

void cell_class::form_conservative(){
    conser[0] = phy.rho;
    conser[1] = phy.rho*phy.u;
    conser[2] = phy.rho*phy.v;
    conser[3] = phy.rho*phy.e;
}

void cell_class::face_normal_out(){
    for(int i=0;i<ecnt;i++){
        fnorm[i] = dot(faces[i]->nor,
                       vec2{faces[i]->mid.x - center.x,faces[i]->mid.y - center.y}) > 0;
    }
}

void cell_class::reform(){
    phy.rho = conser[0];
    phy.u = conser[1]/conser[0];
    phy.v = conser[2]/conser[0];
    phy.e = conser[3]/conser[0];
    phy.T = (phy.e - 0.5*(phy.u*phy.u + phy.v*phy.v))/cc::Cv;
}

void cell_class::copyconver(){
    for(int i=0;i<4;i++){
        conserformer[i] = conser[i];
    }
    tur.miubl_former = tur.miubl;
}

cell_class& gotocell(int number){ return CellList[number-1]; }

face_class& gotoface(int number){ return FaceList[number-1]; }

face_class* link_face(int number){ return &FaceList[number-1]; }

cell_class* link_cell(int number){ return number <= 0 ? nullptr : &CellList[number-1]; }

cell_class* boundary_findcell(face_class* face){
    return face->nei[0] ? face->nei[0] : face->nei[1];
}

bool field_ok(const cell_class& cell){
    const physics& phy = cell.phy;
    return std::isfinite(phy.rho) && phy.rho > 0.0 &&
           std::isfinite(phy.u) && std::isfinite(phy.v) &&
           std::isfinite(phy.T) && phy.T > 0.0 &&
           std::isfinite(phy.a) && std::isfinite(phy.p) && std::isfinite(phy.e) &&
           std::isfinite(cell.tur.miubl) && cell.tur.miubl >= 0.0;
}

void field_mark(const cell_class& cell){
    if(std::isfinite(cell.conser[0]) && cell.conser[0] > 0.0 &&
       std::isfinite(cell.conser[1]) && std::isfinite(cell.conser[2]) &&
       std::isfinite(cell.conser[3]) && std::isfinite(cell.tur.miubl)){
        return;
    }
    int current = field_bad_cell.load(std::memory_order_relaxed);
    while((current == 0 || cell.index < current) &&
          !field_bad_cell.compare_exchange_weak(current,cell.index,
                                                std::memory_order_relaxed)){
    }
}

}
