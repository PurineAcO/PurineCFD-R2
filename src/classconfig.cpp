#include "classconfig.h"
#include "config.h"
#include <cstring>
#include <iostream>

// node_class

namespace cc {

node_class::node_class(int number_,double x_,double y_):
    number(number_),x(x_),y(y_){}

face_class::face_class(int index_,int p1_,int p2_,int c1_,int c2_,short type_):
    index(index_),type(type_),node({p1_,p2_}),cell_1(c1_),cell_2(c2_){
        mid = {0.5 * (NodeList[p1_].x + NodeList[p2_].x),
               0.5 * (NodeList[p1_].y + NodeList[p2_].y)};
        nor = {NodeList[p1_].y - NodeList[p2_].y,
               NodeList[p2_].x - NodeList[p1_].x};
    }

cell_class::cell_class(int index_,int f1_,int f2_,int f3_,int f4_,int ecnt_):
    index(index_),ecnt(ecnt_){
        if(ecnt == 3){face[0] = f1_;face[1] = f2_;face[2] = f3_;face[3] = 0;}
        else if(ecnt == 4){face[0] = f1_;face[1] = f2_;face[2] = f3_;face[3] = f4_;}
        else {std::cerr << "Not Supported Cell Type in" << index << std::endl;}
        if(strncmp(cc::turbulence,"SA", 2) == 0){conser.reserve(5);} 
    }
    
void cell_class::form_conservative(){
    conser = {phy.rho,phy.rho*phy.u,phy.rho*phy.v,phy.rho*phy.e,phy.rho*tur.miubl};
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
    else {return 697;}
}

face_class* link_face(int number){return &gotoface(number);}

cell_class* link_cell(int number){
    if(number == 0){return nullptr;} // 有些边是边界不邻接网格,此时暂时返回空指针.
    return &gotocell(number);}

}