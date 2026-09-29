#include "classconfig.hpp"
#include "config.hpp"
#include "io.hpp"
#include "readmesh.hpp"
#include <cstdio>
#include <iostream>

void print_geom(){
    printf("======GEOM======");
    allcell{
        printf("%d ",icell(i).index);
        // for(int j=0;j<icell(i).ecnt;j++) printf("%d ",icell(i).faces[j]->index);
        printf("\n");
    }
}

int main(){
    // cc::testpath = "debug/debugresult/test.txt";
    // std::freopen(cc::testpath.c_str(),"a", stdout);
    // 加载配置文件
    if(!config::load("config.json"))return 1;
    // 检查网格邻接关系
    if(!readmesh(cc::meshpath.c_str()))return 1;
    print_geom();
}