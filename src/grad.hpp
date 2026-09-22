#pragma once

#include "classconfig.hpp"
#include <cmath>
#include <cstring>
#include <string>

// GGCB梯度
void green_gauss_cell_based(cc::cell_class& cell);
// 面上梯度: 内部面取平均, 壁面修正法向分量(壁面的处理正在进行验证)
void face_gradient(cc::face_class& face);
// 最小二乘LSCB梯度
void least_square_cell_based(cc::cell_class& cell);

static cc::vec2 wall_correct(cc::vec2 grad,double wall,double inner,double nx,double ny,double d){
    double correction = (wall - inner)/d - grad.x*nx - grad.y*ny;
    return cc::vec2{grad.x + correction*nx,grad.y + correction*ny};
}

inline void green_gauss_cell_based(cc::cell_class& cell){
    cell.phy.ugrad = {0.0,0.0};
    cell.phy.vgrad = {0.0,0.0};
    cell.phy.Tgrad = {0.0,0.0};
    cell.tur.miublgrad = {0.0,0.0};
    for(int i=0;i<cell.ecnt;i++){
        const double s = (2*cell.fnorm[i] - 1)*cell.invvol;
        cell.phy.ugrad += (cell.faces[i]->phy.u*s)*cell.faces[i]->nor;
        cell.phy.vgrad += (cell.faces[i]->phy.v*s)*cell.faces[i]->nor;
        cell.phy.Tgrad += (cell.faces[i]->phy.T*s)*cell.faces[i]->nor;
        cell.tur.miublgrad += (cell.faces[i]->tur.miubl*s)*cell.faces[i]->nor;
    }
}

inline void face_gradient(cc::face_class& face){
    if(face.type == cc::INTER){
        face.phy.ugrad = 0.5*(face.nei[0]->phy.ugrad + face.nei[1]->phy.ugrad);
        face.phy.vgrad = 0.5*(face.nei[0]->phy.vgrad + face.nei[1]->phy.vgrad);
        face.phy.Tgrad = 0.5*(face.nei[0]->phy.Tgrad + face.nei[1]->phy.Tgrad);
        face.tur.miublgrad = 0.5*(face.nei[0]->tur.miublgrad + face.nei[1]->tur.miublgrad);
        return;
    }
    face.phy.ugrad = face.phy.vgrad = face.phy.Tgrad = face.tur.miublgrad = {0.0,0.0};
    if(face.type != cc::WALL){
        return;
    }
    cc::cell_class& cell = *cc::boundary_findcell(&face);
    double nx = face.nor.x,ny = face.nor.y;
    double length = std::hypot(nx,ny);
    nx /= length;
    ny /= length;
    double dx = face.mid.x - cell.center.x,dy = face.mid.y - cell.center.y;
    if(nx*dx + ny*dy < 0.0){
        nx = -nx;
        ny = -ny;
    }
    double d = nx*dx + ny*dy;
    face.phy.ugrad = wall_correct(cell.phy.ugrad,0.0,cell.phy.u,nx,ny,d);
    face.phy.vgrad = wall_correct(cell.phy.vgrad,0.0,cell.phy.v,nx,ny,d);
    face.phy.Tgrad = wall_correct(cell.phy.Tgrad,cell.phy.T,cell.phy.T,nx,ny,d);
    face.tur.miublgrad = wall_correct(cell.tur.miublgrad,0.0,cell.tur.miubl,nx,ny,d);
}

inline void least_square_cell_based(cc::cell_class &cell){
    double LU = 0;double DC = 0;double RD = 0;double BU[4] = {};double BD[4] = {};
    for(int i=0;i<cell.ecnt;i++){
        double weight_i,dx_i,dy_i;
        cc::cell_class *nei = cell.nei[i];
        dx_i = nei->center.x - cell.center.x;
        dy_i = nei->center.y - cell.center.y;
        double dis = dx_i * dx_i + dy_i * dy_i;
        weight_i = 1.0/(dis*dis);
        LU += weight_i * dx_i * dx_i;
        DC += weight_i * dx_i * dy_i;
        RD += weight_i * dy_i * dy_i;
        BU[0] += weight_i * dx_i * (nei->phy.u - cell.phy.u);
        BU[1] += weight_i * dx_i * (nei->phy.v - cell.phy.v);
        BU[2] += weight_i * dx_i * (nei->phy.T - cell.phy.T);
        BU[3] += weight_i * dx_i * (nei->tur.miubl - cell.tur.miubl);
        BD[0] += weight_i * dy_i * (nei->phy.u - cell.phy.u);
        BD[1] += weight_i * dy_i * (nei->phy.v - cell.phy.v);
        BD[2] += weight_i * dy_i * (nei->phy.T - cell.phy.T);
        BD[3] += weight_i * dy_i * (nei->tur.miubl - cell.tur.miubl);
    }
    double Det = LU*RD-DC*DC;
    cell.phy.ugrad.x = (RD * BU[0] - DC * BD[0])/Det;
    cell.phy.ugrad.y = (LU * BD[0] - DC * BU[0])/Det;
    cell.phy.vgrad.x = (RD * BU[1] - DC * BD[1])/Det;
    cell.phy.vgrad.y = (LU * BD[1] - DC * BU[1])/Det;
    cell.phy.Tgrad.x = (RD * BU[2] - DC * BD[2])/Det;
    cell.phy.Tgrad.y = (LU * BD[2] - DC * BU[2])/Det;
    cell.tur.miublgrad.x = (RD * BU[3] - DC * BD[3])/Det;
    cell.tur.miublgrad.y = (LU * BD[3] - DC * BU[3])/Det;
}