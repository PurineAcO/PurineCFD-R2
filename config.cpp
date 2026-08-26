#include <cstdlib>
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

}