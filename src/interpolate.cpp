#include "interpolat.h"
#include "classconfig.h"

void interpolate_mid(cc::cell_class& cell){
    for(int i=0;i<cell.ecnt;i++){
        cc::face_class& face = *cell.faces[i];
        if(face.nei[0] == &cell){
            face.face_physic_mid();
        }
    }
}
