#pragma once

#include "classconfig.hpp"
#include "config.hpp"
#include <cmath>
#include <linux/stat.h>

// 无粘通量
void convect_JST(cc::face_class& face);
// 汇总单元各面的无粘和黏性通量
void assemble_flux(cc::cell_class& cell);
// Roe无粘通量
void convect_ROE(cc::face_class &face);

inline void convect_JST(cc::face_class& face){
    const double rho = face.phy.rho,u = face.phy.u,v = face.phy.v;
    const double H = cc::Cp*face.phy.T + 0.5*(u*u + v*v);
    const double p = cc::R*rho*face.phy.T;
    const double F[4] = {rho*u,rho*u*u + p,rho*u*v,rho*u*H};
    const double G[4] = {rho*v,rho*u*v,rho*v*v + p,rho*v*H};
    for(int j=0;j<4;j++){
        face.convect[j] = F[j]*face.nor.x + G[j]*face.nor.y;
    }
}

inline void assemble_flux(cc::cell_class& cell){
    for(int j=0;j<4;j++){
        cell.convect[j] = 0.0;
        cell.visflux[j] = 0.0;
    }
    for(int i=0;i<cell.ecnt;i++){
        cc::face_class& face = *cell.faces[i];
        const int outer = 2*cell.fnorm[i] - 1;
        for(int j=0;j<4;j++){
            cell.convect[j] += outer*face.convect[j];
            cell.visflux[j] += outer*face.visflux[j];
        }
    }
}

inline void convect_ROE(cc::face_class &face){
    cc::cell_class *high = face.high;cc::cell_class *low = face.low;
    // Roe平均值
    double Hl = cc::Cp * low->phy.T + 0.5 * (low->phy.u * low->phy.u + low->phy.v * low->phy.v);
    double Hh = cc::Cp * high->phy.T + 0.5 * (high->phy.u * high->phy.u + high->phy.v * high->phy.v);
    double rhobl = std::sqrt(low->phy.rho * high->phy.rho);
    double ubl = (std::sqrt(low->phy.rho)*low->phy.u + std::sqrt(high->phy.rho)*high->phy.u)/
                  (std::sqrt(low->phy.rho)+std::sqrt(high->phy.rho));
    double vbl = (std::sqrt(low->phy.rho)*low->phy.v + std::sqrt(high->phy.rho)*high->phy.v)/
                  (std::sqrt(low->phy.rho)+std::sqrt(high->phy.rho));
    double Hbl = (std::sqrt(low->phy.rho)*Hl + std::sqrt(high->phy.rho)*Hh)/
                  (std::sqrt(low->phy.rho)+std::sqrt(high->phy.rho));
    double abl = std::sqrt((cc::gamma-1) * (Hbl - 0.5 * (ubl*ubl + vbl*vbl)));
    double length = std::sqrt(face.nor.x*face.nor.x + face.nor.y*face.nor.y);
    double nx = face.nor.x/length;double ny = face.nor.y/length;
    double unbl = ubl * nx + vbl * ny;
    double utbl = ubl * ny - vbl * nx;

    // 特征值
    double lambda[4] = {unbl - abl,unbl,unbl,unbl + abl}; 

    // 波强度
    double dp = cc::R * (high->phy.T * high->phy.rho - low->phy.T * low->phy.rho);
    double drho = high->phy.rho - low->phy.rho;
    double dun = nx * (high->phy.u - low->phy.u) + ny * (high->phy.v - low->phy.v);
    double dut = ny * (high->phy.u - low->phy.u) - nx * (high->phy.v - low->phy.v);
    double alpha[4] = {(dp-rhobl*abl*dun)/(2*abl*abl),drho-dp/(2*abl*abl),rhobl*dut,(dp+rhobl*abl*dun)/(2*abl*abl)};
    
    // 特征向量
    double tz[4][4] = {{1,ubl-abl*nx,vbl-abl*ny,Hbl-unbl*abl},
                       {1,ubl,       vbl,       (ubl*ubl+vbl*vbl)*0.5},
                       {0,-abl*ny,   abl*nx,    -abl*utbl},
                       {1,ubl+abl*nx,vbl+abl*ny,Hbl+unbl*abl}};

    // 左右无粘通量
    double unl = low->phy.u * nx + low->phy.v * ny;
    double unh = high->phy.u * nx + high->phy.v * ny;
    double FL[4] = {low->phy.rho * unl ,low->phy.rho*low->phy.u*unl+low->phy.p*nx,
                    low->phy.rho*low->phy.v*unl+low->phy.p*ny,unl*(low->phy.rho*low->phy.e+low->phy.p)};
    double FR[4] = {high->phy.rho * unh ,high->phy.rho*high->phy.u*unh+high->phy.p*nx,
                    high->phy.rho*high->phy.v*unh+high->phy.p*ny,unh*(high->phy.rho*high->phy.e+high->phy.p)};

    // 形成面上对流通量
    vecfor(4) face.convect[i] = 0;
    vecfor(4) face.convect[i] += 0.5*(FL[i]+FR[i]);
    vecfor(4) for(int j=0;j<4;j++) face.convect[j] -= 0.5 * lambda[i] * alpha[i] * tz[i][j];
        
}