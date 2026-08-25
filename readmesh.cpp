#include <cstdio>
#include <cstring>
#include <iostream>
#include <ostream>
#include "classconfig.h"

int main(int argc,char* argv[]){

    const char* path = (argc > 1) ? argv[1] : "mesh/tunnel.txt";
    FILE* fp = fopen(path, "rb");
    if(fp == NULL){std::cerr << "Cannot open the mesh" << std::endl;return 1;}

    int cell_num,edge_num,node_num;
    if(fscanf(fp,"%d %d %d",&node_num,&edge_num,&cell_num) != 3){
        std::cerr << "Not Correct Mesh" << std::endl;return 1;
    }
    std::cout << "Node Nums:" << node_num << std::endl;
    std::cout << "Edge Nums:" << edge_num << std::endl;
    std::cout << "Cell Nums:" << cell_num << std::endl;

    char line[100];fgets(line,sizeof(line),fp); // 吃掉一个换行符
    
    if(fgets(line,sizeof(line),fp) == NULL){
        std::cerr << "No Node Line" << std::endl;return 1;
    }
    if(strncmp(line,"(node)",6) != 0){
        std::cerr << "Bad Format for Node" << std::endl;return 1;
    }

    for(int i=1;i<=node_num;i++){
        int number;double x,y;
        if(fscanf(fp,"%d %lf %lf",&number,&x,&y) != 3){
            std::cerr << "Node data read error at line " << i << std::endl;return 1;
        }
        if(i!=number){std::cerr << "Skip Node" << std::endl;}
        cc::node_class(number,x,y);
    }

    printf("Node Read OK\n");

    fclose(fp);

    return 0;
    
}