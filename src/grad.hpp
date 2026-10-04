#pragma once

#include "classconfig.hpp"
#include <cmath>
#include <cstring>

// GGCB梯度
void green_gauss_cell_based(cc::cell_class& cell);
// LSCB梯度
void least_square_cell_based(cc::cell_class& cell);

inline void green_gauss_cell_based(cc::cell_class& cell) {
    cell.phgrad.clear();
    cell.tur.miublgrad = {0.0, 0.0};
    allface(cell) {
        const double s = (2 * cell.fnorm[i] - 1) / cell.vol;
        cell.phgrad = cell.phgrad + (cell.faces[i]->phy * s) * cell.faces[i]->nor;
        cell.tur.miublgrad += (cell.faces[i]->tur.miubl * s) * cell.faces[i]->nor;
    }
}

inline void least_square_cell_based(cc::cell_class& cell) {
    // 预处理过程在least_square_cell_based_preprocess里,该函数定义在geometry
    double BU[5] = {};
    double BD[5] = {};
    allface(cell) {
        cc::cell_class* nei = cell.nei[i];
        // 虚网格不带几何, 不参与最小二乘
        if (nei == nullptr || nei->index > cc::cell_num) continue;
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
    double Det = cell.LSCB.LU * cell.LSCB.RD - cell.LSCB.DC * cell.LSCB.DC;
    cell.phgrad.ugrad.x = (cell.LSCB.RD * BU[0] - cell.LSCB.DC * BD[0]) / Det;
    cell.phgrad.ugrad.y = (cell.LSCB.LU * BD[0] - cell.LSCB.DC * BU[0]) / Det;
    cell.phgrad.vgrad.x = (cell.LSCB.RD * BU[1] - cell.LSCB.DC * BD[1]) / Det;
    cell.phgrad.vgrad.y = (cell.LSCB.LU * BD[1] - cell.LSCB.DC * BU[1]) / Det;
    cell.phgrad.Tgrad.x = (cell.LSCB.RD * BU[2] - cell.LSCB.DC * BD[2]) / Det;
    cell.phgrad.Tgrad.y = (cell.LSCB.LU * BD[2] - cell.LSCB.DC * BU[2]) / Det;
    cell.tur.miublgrad.x = (cell.LSCB.RD * BU[3] - cell.LSCB.DC * BD[3]) / Det;
    cell.tur.miublgrad.y = (cell.LSCB.LU * BD[3] - cell.LSCB.DC * BU[3]) / Det;
    cell.phgrad.rhograd.x = (cell.LSCB.RD * BU[4] - cell.LSCB.DC * BD[4]) / Det;
    cell.phgrad.rhograd.y = (cell.LSCB.LU * BD[4] - cell.LSCB.DC * BU[4]) / Det;
}