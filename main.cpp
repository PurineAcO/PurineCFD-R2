#include <cstdio>
#include "readmesh.h"
#include "geometry.h"
#include "config.h"

int main(){
    if(readmesh(cc::path)){return 1;}
    if(linkmesh()){return 1;}
    geometry();
    return 0;
}