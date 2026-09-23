#pragma once

#include "classconfig.hpp"
#include "readmesh.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

// 找到全部邻接节点
void findnode(cc::cell_class& cell);
// 计算网格体积(四个三角形面积和的一半)
void volume(cc::cell_class& cell);
// 计算网格质心(按三角形面积加权)
void center(cc::cell_class& cell);
// 计算壁面距离
void sad(cc::cell_class& cell);
// 最小二乘LSCB梯度预处理
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
    for(int i=0;i<cell.ecnt;i++){
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
    for(int i=0;i<cell.ecnt;i++){
        cc::cell_class *nei = cell.nei[i];
        cell.LSCB.dxi[i] = nei->center.x - cell.center.x;
        cell.LSCB.dyi[i] = nei->center.y - cell.center.y;
        cell.LSCB.wi[i] = 1/(cell.LSCB.dxi[i]*cell.LSCB.dxi[i] + cell.LSCB.dyi[i]*cell.LSCB.dyi[i]);
        cell.LSCB.LU += cell.LSCB.wi[i] * cell.LSCB.dxi[i] * cell.LSCB.dxi[i];
        cell.LSCB.DC += cell.LSCB.wi[i] * cell.LSCB.dxi[i] * cell.LSCB.dyi[i];
        cell.LSCB.RD += cell.LSCB.wi[i] * cell.LSCB.dyi[i] * cell.LSCB.dyi[i];
    }
}

inline bool geometrymain(){
    if(!linkmesh()){
        return false;
    }
    // 带结构化邻接表时, 在面邻接已建立后补齐东西南北四个方向
    if(structer::ifstructer && !link_structed_mesh()){
        return false;
    }
    for(cc::cell_class& cell : cc::CellList){
        findnode(cell);
        volume(cell);
        center(cell);
        cell.face_normal_out();
        if(!std::isfinite(cell.vol) || cell.vol <= 0.0){
            fprintf(stderr,"Error: cell #%d has an invalid area\n",cell.index);
            return false;
        }
        cell.invvol = 1.0/cell.vol;
        for(int i=0;i<cell.ecnt;i++){
            cell.proj.x += 0.5*std::abs(cell.faces[i]->nor.x);
            cell.proj.y += 0.5*std::abs(cell.faces[i]->nor.y);
        }
        least_square_cell_based_preprocess(cell);
    }
    return true;
}
