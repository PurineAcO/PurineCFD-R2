#include <cstdio>
#include <cstring>
#include <iostream>
#include <ostream>
#include <string>
#include <vector>
#include <algorithm>
#include "classconfig.h"
#include "config.h"

short readmesh(const char* path){

    FILE* fp = fopen(path, "rb");
    if(fp == NULL){std::cerr << "Cannot open the mesh" << std::endl;return 1;}

    int node_num;
    std::vector<std::string>type_name,type_std;

    if(fscanf(fp,"%d %d %d %d",
        &node_num,&cc::face_num,&cc::cell_num,&cc::type_total) != 4){
        std::cerr << "Not Correct Mesh" << std::endl;return 1;
    }
    printf("Node:%d,Face:%d,Cell:%d,Types:%d\n",
                    node_num,cc::face_num,cc::cell_num,cc::type_total);

    char line[100];fgets(line,sizeof(line),fp); // 吸收换行符
    type_name.reserve(cc::type_total);
    type_std.reserve(cc::type_total);
    for(int i = 1;i<=cc::type_total;i++){
        char type_aka[50],type_std_name[50];
        fscanf(fp,"%49s = %49s",type_aka,type_std_name);
        type_name.push_back(type_aka);
        type_std.push_back(type_std_name);
    }

    fgets(line,sizeof(line),fp); // 吸收换行符

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

    if (cc::NodeList.size() == (size_t)node_num){printf("Node Read OK\n");}
    else{std::cerr << "Node Wrong Nums" << std::endl;}
    
    fgets(line,sizeof(line),fp); // 吸收换行符
    fgets(line,sizeof(line),fp); // 吸收(end node),没有这一行一般不报错.
    if(fgets(line,sizeof(line),fp) == NULL){
        std::cerr << "No Edge Line" << std::endl;return 1;
    }
    if(strncmp(line,"(edge)",6) != 0){
        std::cerr << "Bad Format for Edge" << std::endl;return 1;
    }

    cc::FaceList.reserve(cc::face_num);
    for (int s=1;s<=cc::type_total;s++){
        if(fgets(line,sizeof(line),fp) == NULL){
            std::cerr << "No Face Type Line" << std::endl;return 1;
        }
        line[strcspn(line,"\r\n")] = 0;   // 去掉换行符, 便于 strcmp
        if(strcmp(line,type_name[s-1].c_str())!=0){
            std::cerr << "Not Correct Face Type: got '" << line
                      << "' expect '" << type_name[s-1] << "'" << std::endl;return 1;
        }
        short type_code = cc::get_facetype(type_std[s-1].c_str());
        if(type_code == 697){
            std::cerr << "Not Supported Face Type: '" << type_std[s-1] << "'" << std::endl;return 1;
        }
        while(true){
            if(fgets(line,sizeof(line),fp) == NULL){
                std::cerr << "Edge data read error in section " << s << std::endl;return 1;
            }
            if(strncmp(line,"(end",4) == 0) break;   // 段结束标记
            int index,p1,p2,c1,c2;
            if(sscanf(line,"%d%d%d%d%d",&index,&p1,&p2,&c1,&c2) != 5){
                std::cerr << "Edge data read error in section " << s << std::endl;return 1;
            }
            cc::FaceList.emplace_back(index,p1,p2,c1,c2,type_code);
        }
        printf("Read OK:%s\n",type_name[s-1].c_str());
    }

    std::sort(cc::FaceList.begin(), cc::FaceList.end(),
              [](const cc::face_class& a, const cc::face_class& b)
              { return a.index < b.index; });

    if(cc::FaceList.size() != (size_t)cc::face_num){
        std::cerr << "Edge Wrong Nums" << std::endl;return 1;
    }
    
    printf("Edge Read OK\n");

    // fgets(line,sizeof(line),fp); // 吸收换行符
    fgets(line,sizeof(line),fp); // 吸收(end edge),没有这一行一般不报错.
    if(fgets(line,sizeof(line),fp) == NULL){
        std::cerr << "No Cell Line" << std::endl;return 1;
    }
    if(strncmp(line,"(cell)",6) != 0){
        std::cerr << "Bad Format for Cell" << std::endl;return 1;
    }

    cc::CellList.reserve(cc::cell_num);
    for(int i=1;i<=cc::cell_num;i++){
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

    if(cc::CellList.size() != (size_t)cc::cell_num){
        std::cerr << "Cell Wrong Nums" << std::endl;return 1;
    }
    
    for(int i=0;i<(int)cc::CellList.size();i++){
        if(cc::CellList[i].index != i+1){
            std::cerr << "Cell order broken at position " << i << std::endl;return 1;
        }
    }

    printf("Cell Read OK\n");

    fclose(fp);
    
    return 0;
}

short linkmesh(){

    for(cc::cell_class& cell: cc::CellList){
        for(int i = 0; i<cell.ecnt ; i++){
            cell.nei[i] = cc::link_face(cell.face[i]);
        }
    }

    for(cc::face_class& face : cc::FaceList){
        face.nei[0] = cc::link_cell(face.cell_1);
        face.nei[1] = cc::link_cell(face.cell_2);
        if(face.type == cc::WALL){
            cc::WallFaces.push_back(&face);
        }else if(face.type == cc::VIL){
            cc::VILFaces.push_back(&face);
        }else if(face.type == cc::POL){
            cc::POLFaces.push_back(&face);
        }
    }

    printf("Link Cell and Face OK \n");
    printf("Wall:%zu, Inlet:%zu, Outlet:%zu\n",
           cc::WallFaces.size(), cc::VILFaces.size(), cc::POLFaces.size());
    return 0;
}