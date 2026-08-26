#include <cstdio>
#include "readmesh.h"
#include "config.h"

int main(){
    if(readmesh(cc::path)){return 1;}
    return 0;
}