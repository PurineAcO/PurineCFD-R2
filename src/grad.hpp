#pragma once

#include "classconfig.hpp"
#include <cmath>
#include <cstring>

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
    // 预处理过程在least_square_cell_based_preprocess里,该函数定义在geometry
    double BU[5] = {};double BD[5] = {};
    for(int i=0;i<cell.ecnt;i++){
        cc::cell_class *nei = cell.nei[i];
        if(nei == nullptr)continue;
        BU[0] += cell.LSCB.wi[i] * cell.LSCB.dxi[i] * (nei->phy.u - cell.phy.u);
        BU[1] += cell.LSCB.wi[i] * cell.LSCB.dxi[i] * (nei->phy.v - cell.phy.v);
        BU[2] += cell.LSCB.wi[i] * cell.LSCB.dxi[i] * (nei->phy.T - cell.phy.T);
        BU[3] += cell.LSCB.wi[i] * cell.LSCB.dxi[i] * (nei->tur.miubl - cell.tur.miubl);
        BU[4] += cell.LSCB.wi[i] * cell.LSCB.dxi[i] * (nei->phy.rho - cell.phy.rho);
        BD[0] += cell.LSCB.wi[i] * cell.LSCB.dyi[i] * (nei->phy.u - cell.phy.u);
        BD[1] += cell.LSCB.wi[i] * cell.LSCB.dyi[i] * (nei->phy.v - cell.phy.v);
        BD[2] += cell.LSCB.wi[i] * cell.LSCB.dyi[i] * (nei->phy.T - cell.phy.T);
        BD[3] += cell.LSCB.wi[i] * cell.LSCB.dyi[i] * (nei->tur.miubl - cell.tur.miubl);
        BD[4] += cell.LSCB.wi[i] * cell.LSCB.dyi[i] * (nei->phy.rho - cell.phy.rho);
    }
    double Det = cell.LSCB.LU*cell.LSCB.RD-cell.LSCB.DC*cell.LSCB.DC;
    cell.phy.ugrad.x = (cell.LSCB.RD * BU[0] - cell.LSCB.DC * BD[0])/Det;
    cell.phy.ugrad.y = (cell.LSCB.LU * BD[0] - cell.LSCB.DC * BU[0])/Det;
    cell.phy.vgrad.x = (cell.LSCB.RD * BU[1] - cell.LSCB.DC * BD[1])/Det;
    cell.phy.vgrad.y = (cell.LSCB.LU * BD[1] - cell.LSCB.DC * BU[1])/Det;
    cell.phy.Tgrad.x = (cell.LSCB.RD * BU[2] - cell.LSCB.DC * BD[2])/Det;
    cell.phy.Tgrad.y = (cell.LSCB.LU * BD[2] - cell.LSCB.DC * BU[2])/Det;
    cell.tur.miublgrad.x = (cell.LSCB.RD * BU[3] - cell.LSCB.DC * BD[3])/Det;
    cell.tur.miublgrad.y = (cell.LSCB.LU * BD[3] - cell.LSCB.DC * BU[3])/Det;
    cell.phy.rhograd.x = (cell.LSCB.RD * BU[4] - cell.LSCB.DC * BD[4])/Det;
    cell.phy.rhograd.y = (cell.LSCB.LU * BD[4] - cell.LSCB.DC * BU[4])/Det;
}