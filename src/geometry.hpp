#pragma once

#include "classconfig.hpp"
#include "readmesh.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

// 找到全部邻接节点
void findnode(cc::cell_class& cell);
// 计算网格体积
void volume(cc::cell_class& cell);
// 计算网格质心
void center(cc::cell_class& cell);
// 计算壁面距离
void sad(cc::cell_class& cell);
// LSCB梯度预处理
void least_square_cell_based_preprocess(cc::cell_class &cell);
// 几何分析主程序
bool geometrymain();

static bool contains_node(int num[],int number){
    for(int i=0;i<4;i++){
        if(number == num[i]){
            return true;
        }
    }
    return false;
}

static double point_distance(cc::vec2 a,cc::vec2 b){
    return std::sqrt((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y));
}

static double triangle_area(double x1,double y1,double x2,double y2,double x3,double y3){
    return std::abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))*0.5;
}

static cc::vec2 triangle_center(double x1,double y1,double x2,double y2,double x3,double y3){
    return {(x1+x2+x3)/3,(y1+y2+y3)/3};
}

inline void findnode(cc::cell_class& cell){
    short place = 0;
    allface(cell){
        if(!contains_node(cell.node,cell.faces[i]->node[0]->number - 1)){
            cell.node[place] = cell.faces[i]->node[0]->number - 1;
            place++;
        }
        if(!contains_node(cell.node,cell.faces[i]->node[1]->number - 1)){
            cell.node[place] = cell.faces[i]->node[1]->number - 1;
            place++;
        }
        if(place == cell.ecnt){
            break;
        }
    }
}

inline void volume(cc::cell_class& cell){
    const cc::node_class& p0 = cc::NodeList[cell.node[0]];
    const cc::node_class& p1 = cc::NodeList[cell.node[1]];
    const cc::node_class& p2 = cc::NodeList[cell.node[2]];
    const cc::node_class& p3 = cc::NodeList[cell.node[3]];
    cell.vol = (triangle_area(p0.x,p0.y,p1.x,p1.y,p2.x,p2.y) +
                triangle_area(p0.x,p0.y,p1.x,p1.y,p3.x,p3.y) +
                triangle_area(p0.x,p0.y,p2.x,p2.y,p3.x,p3.y) +
                triangle_area(p1.x,p1.y,p2.x,p2.y,p3.x,p3.y))*0.5;
}

inline void center(cc::cell_class& cell){
    const cc::node_class& p0 = cc::NodeList[cell.node[0]];
    const cc::node_class& p1 = cc::NodeList[cell.node[1]];
    const cc::node_class& p2 = cc::NodeList[cell.node[2]];
    const cc::node_class& p3 = cc::NodeList[cell.node[3]];
    double S[4];
    cc::vec2 G[4];
    S[0] = triangle_area(p0.x,p0.y,p1.x,p1.y,p2.x,p2.y);
    G[0] = triangle_center(p0.x,p0.y,p1.x,p1.y,p2.x,p2.y);
    S[1] = triangle_area(p0.x,p0.y,p1.x,p1.y,p3.x,p3.y);
    G[1] = triangle_center(p0.x,p0.y,p1.x,p1.y,p3.x,p3.y);
    S[2] = triangle_area(p0.x,p0.y,p2.x,p2.y,p3.x,p3.y);
    G[2] = triangle_center(p0.x,p0.y,p2.x,p2.y,p3.x,p3.y);
    S[3] = triangle_area(p1.x,p1.y,p2.x,p2.y,p3.x,p3.y);
    G[3] = triangle_center(p1.x,p1.y,p2.x,p2.y,p3.x,p3.y);
    cell.center = {(S[0]*G[0].x + S[1]*G[1].x + S[2]*G[2].x + S[3]*G[3].x)/
                       (S[0] + S[1] + S[2] + S[3]),
                   (S[0]*G[0].y + S[1]*G[1].y + S[2]*G[2].y + S[3]*G[3].y)/
                       (S[0] + S[1] + S[2] + S[3])};
}

inline void sad(cc::cell_class& cell){
    double distance = std::numeric_limits<double>::infinity();
    for(cc::face_class* wall : cc::WallFaces){
        double length = point_distance(cell.center,wall->mid);
        distance = std::min(distance,length);
    }
    cell.tur.sad = distance;
}

inline void least_square_cell_based_preprocess(cc::cell_class &cell){
    allface(cell){
        cc::cell_class *nei = cell.nei[i];
        // 虚网格不带几何, 不参与最小二乘
        if(nei == nullptr || nei->index > cc::cell_num){
            continue;
        }
        cell.LSCB.dxi[i] = nei->center.x - cell.center.x;
        cell.LSCB.dyi[i] = nei->center.y - cell.center.y;
        cell.LSCB.wi[i] = 1/(cell.LSCB.dxi[i]*cell.LSCB.dxi[i] + cell.LSCB.dyi[i]*cell.LSCB.dyi[i]);
        cell.LSCB.LU += cell.LSCB.wi[i] * cell.LSCB.dxi[i] * cell.LSCB.dxi[i];
        cell.LSCB.DC += cell.LSCB.wi[i] * cell.LSCB.dxi[i] * cell.LSCB.dyi[i];
        cell.LSCB.RD += cell.LSCB.wi[i] * cell.LSCB.dyi[i] * cell.LSCB.dyi[i];
    }
}

inline bool geometrymain(){
    if(!linkmesh())return false;
    if(structer::ifstructer && !link_structed_mesh())return false;
    allcell{
        findnode(icell(i));
        volume(icell(i));
        center(icell(i));
        icell(i).face_normal_out();
        sad(icell(i));
    }
    // face.outer 依赖格子质心, 必须在 center 之后才定
    for(cc::face_class& face : cc::FaceList){
        if(face.type == cc::INTER) face.normal_out();
    }
    return true;
}
