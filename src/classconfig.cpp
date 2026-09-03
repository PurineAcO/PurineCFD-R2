#include "classconfig.h"
#include "config.h"
#include "physic.h"
#include <cmath>
#include <iostream>

namespace cc {

node_class::node_class(int number_,double x_,double y_):
    number(number_),x(x_),y(y_){}

face_class::face_class(int index_,int p1_,int p2_,int c1_,int c2_,short type_):
    index(index_),type(type_),node({p1_,p2_}),cell_1(c1_),cell_2(c2_){
        mid = {0.5 * (NodeList[p1_-1].x + NodeList[p2_-1].x),
               0.5 * (NodeList[p1_-1].y + NodeList[p2_-1].y)};
        nor = {NodeList[p1_-1].y - NodeList[p2_-1].y,
               NodeList[p2_-1].x - NodeList[p1_-1].x};
    }

cell_class::cell_class(int index_,int f1_,int f2_,int f3_,int f4_,int ecnt_):
    index(index_),ecnt(ecnt_){
        if(ecnt == 3){face[0] = f1_;face[1] = f2_;face[2] = f3_;face[3] = 0;}
        else if(ecnt == 4){face[0] = f1_;face[1] = f2_;face[2] = f3_;face[3] = f4_;}
        else {std::cerr << "Not Supported Cell Type in" << index << std::endl;}
    }

void cell_class::form_conservative(){
    conser[0] = phy.rho;
    conser[1] = phy.rho*phy.u;
    conser[2] = phy.rho*phy.v;
    conser[3] = phy.rho*phy.e;
}

void cell_class::face_normal_out(){
    for(int i=0;i<ecnt;i++){
        fnorm[i] = dot(nei[i]->nor, {nei[i]->mid.x - center.x, nei[i]->mid.y - center.y}) > 0;
    }
}

void cell_class::form_physic(){
    phy.e = cc::Cv * phy.T + 0.5*(phy.u*phy.u+phy.v*phy.v);
    phy.p = cc::R * phy.rho * phy.T;
    phy.a = get_sonic_velocity(phy.T);
}

void cell_class::reform(){
    phy.rho = conser[0];
    phy.u = conser[1]/conser[0];
    phy.v = conser[2]/conser[0];
    phy.e = conser[3]/conser[0];
    phy.T = (phy.e - 0.5*(phy.u*phy.u+phy.v*phy.v))/cc::Cv;
}

void cell_class::copyconver(){for(int i=0;i<NEQ;i++){conserformer[i] = conser[i];}}

void face_class::face_physic_mid(){
    if(type == INTER){
        phy.rho = 0.5 * (nei[0]->phy.rho + nei[1]->phy.rho);
        phy.u = 0.5 * (nei[0]->phy.u + nei[1]->phy.u);
        phy.v = 0.5 * (nei[0]->phy.v + nei[1]->phy.v);
        phy.T = 0.5 * (nei[0]->phy.T + nei[1]->phy.T);
    }
}

void face_class::form_physic(){
    phy.e = cc::Cv * phy.T + 0.5*(phy.u*phy.u+phy.v*phy.v);
    phy.p = cc::R * phy.rho * phy.T;
    phy.a = get_sonic_velocity(phy.T);
}

cell_class& gotocell(int number){
    if(number < 1 || number > cc::cell_num){
        std::cerr << "Cell Not Exists" << std::endl;
        std::exit(1);}
    return cc::CellList[number - 1];
}

face_class& gotoface(int number){
    if(number < 1 || number > cc::face_num){
        std::cerr << "Edge Not Exists" << std::endl;
        std::exit(1);
    }
    return cc::FaceList[number - 1];
}

short get_facetype(const char *face_std_name){
    if(strcmp(face_std_name,"WALL") == 0){return WALL;}
    else if(strcmp(face_std_name,"INTER") == 0){return INTER;}
    else if(strcmp(face_std_name,"VIL") == 0){return VIL;}
    else if(strcmp(face_std_name,"POL") == 0){return POL;}
    else if(strcmp(face_std_name,"FAR") == 0){return FAR;}
    else {return 697;}
}

face_class* link_face(int number){return &gotoface(number);}

cell_class* link_cell(int number){
    if(number == 0){return nullptr;}
    return &gotocell(number);
}

cell_class* find_neicell(int index,face_class* face){
    if(face->nei[0] == NULL || face->nei[1] == NULL){return nullptr;}
    return (face->nei[0]->index == index) ? face->nei[1] : face->nei[0];
}

}
