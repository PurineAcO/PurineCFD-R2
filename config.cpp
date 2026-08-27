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

}