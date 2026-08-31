#include "classconfig.h"

void interpolate_mid(cc::cell_class &cell){
    for(int i=0;i<cell.ecnt;i++){
        cell.nei[i]->face_physic_mid();
        cell.nei[i]->form_physic();
    }
}