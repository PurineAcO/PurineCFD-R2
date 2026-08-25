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

    char line[100];fgets(line,sizeof(line),fp); // 吸收换行符
    
    if(fgets(line,sizeof(line),fp) == NULL){
        std::cerr << "No Node Line" << std::endl;return 1;
    }
    if(strncmp(line,"(node)",6) != 0){
        std::cerr << "Bad Format for Node" << std::endl;return 1;
    }

    cc::NodeList.reserve(node_num);
    for(int i=1;i<=node_num;i++){
        int number;double x,y;
        if(fscanf(fp,"%d %lf %lf",&number,&x,&y) != 3){
            std::cerr << "Node data read error at line " << i << std::endl;return 1;
        }
        if(i!=number){std::cerr << "Skip Node" << std::endl;}
        cc::NodeList.emplace_back(number,x,y);
    }

    printf("Node Read OK\n");

    fgets(line,sizeof(line),fp); // 吸收换行符
    fgets(line,sizeof(line),fp); // 吸收(end node),没有这一行一般不报错.
    if(fgets(line,sizeof(line),fp) == NULL){
        std::cerr << "No Edge Line" << std::endl;return 1;
    }
    if(strncmp(line,"(edge)",6)){
        std::cerr << "Bad Format for Edge" << std::endl;
    }

    cc::FaceList.reserve(edge_num);
    for(int i=1;i<=edge_num;i++){
        int index,p1,p2,c1,c2;
        if(fscanf(fp,"%d%d%d%d%d",&index,&p1,&p2,&c1,&c2) != 5){
            std::cerr << "Edge data read error at line " << i << std::endl;return 1;
        }
        if(i!=index){std::cerr << "Skip Edge" << std::endl;}
        cc::FaceList.emplace_back(index,p1,p2,c1,c2);
    }

    printf("Edge Read OK\n");

    fgets(line,sizeof(line),fp); // 吸收换行符
    fgets(line,sizeof(line),fp); // 吸收(end edge),没有这一行一般不报错.
    if(fgets(line,sizeof(line),fp) == NULL){
        std::cerr << "No Cell Line" << std::endl;return 1;
    }
    if(strncmp(line,"(cell)",6)){
        std::cerr << "Bad Format for Cell" << std::endl;
    }

    cc::CellList.reserve(cell_num);
    for(int i=1;i<=cell_num;i++){
        if(fgets(line,sizeof(line),fp) == NULL){
            std::cerr << "Cell data read error at line " << i << std::endl;return 1;
        }
        int index, f1=0, f2=0, f3=0, f4=0;
        int n = sscanf(line, "%d %d %d %d %d", &index,&f1,&f2,&f3,&f4);
        if((n-1) < 2 || (n-1) > 4){
            std::cerr << "Bad cell line:Incorrect Nums of Edge at" << i << std::endl;return 1;
        }
        cc::CellList.emplace_back(index, f1, f2, f3, f4, n-1);
    }

    printf("Cell Read OK\n");
    fclose(fp);
    
    return 0;
    
}