#include <cstdlib>
#include <cstring>
#include <iostream>
#include "config.h"
#include "classconfig.h"

namespace cc {

cc::cell_class& gotocell(int number){
    if(number < 1 || number > cc::cell_num){
        std::cerr << "Cell Not Exists" << std::endl;
        std::exit(1);}
    return cc::CellList[number - 1];
}

cc::face_class& gotoface(int number){
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

cc::face_class* link_face(int number){return &gotoface(number);}

cc::cell_class* link_cell(int number){
    if(number == 0){return nullptr;} // 有些边是边界不邻接网格,此时暂时返回空指针.
    return &gotocell(number);}

}

// 用户自定义部分

void DO_BEFORE_SOLVE(){
    // 边界条件
    cc::Ma = 0.3;cc::T = 300;cc::Re = 3e6;cc::AOA = 0;cc::p = 1e5;
}