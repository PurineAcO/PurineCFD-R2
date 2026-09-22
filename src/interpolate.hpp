#pragma once

#include "classconfig.hpp"

// 面上中心差分插值; 内部面只由nei[0]侧单元调用
void interpolate_mid(cc::cell_class& cell);

inline void interpolate_mid(cc::cell_class& cell){
    for(int i=0;i<cell.ecnt;i++){
        cc::face_class& face = *cell.faces[i];
        if(face.nei[0] == &cell){
            face.face_physic_mid();
        }
    }
}
