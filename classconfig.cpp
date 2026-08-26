#include "classconfig.h"
#include <iostream>

// node_class

namespace cc {

node_class::node_class(int number_,double x_,double y_):
    number(number_),x(x_),y(y_){}

face_class::face_class(int index_,int p1_,int p2_,int c1_,int c2_):
    index(index_),cell_1(c1_),cell_2(c2_){
        mid = {0.5 * (NodeList[p1_].x + NodeList[p2_].x),
               0.5 * (NodeList[p1_].y + NodeList[p2_].y)};
        nor = {NodeList[p1_].y - NodeList[p2_].y,
               NodeList[p2_].x - NodeList[p1_].x};
    }

cell_class::cell_class(int index_,int f1_,int f2_,int f3_,int f4_,int ecnt_):
    index(index_),ecnt(ecnt_){
        if(ecnt == 3){face[0] = f1_;face[1] = f2_;face[2] = f3_;face[3] = 0;}
        else if(ecnt == 4){face[0] = f1_;face[2] = f2_;face[3] = f3_;face[3] = f4_;}
        else {std::cerr << "Not Supported Cell Type in" << index << std::endl;}
    }


}