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
    // MUSCL重构后的左右状态
    const cc::physics& L = face.lowp;
    const cc::physics& R = face.highp;
    // Roe平均值
    const double sl = std::sqrt(L.rho),sh = std::sqrt(R.rho);
    const double Hl = cc::Cp * L.T + 0.5 * (L.u * L.u + L.v * L.v);
    const double Hh = cc::Cp * R.T + 0.5 * (R.u * R.u + R.v * R.v);
    const double rhobl = sl*sh;
    const double ubl = (sl*L.u + sh*R.u)/(sl+sh);
    const double vbl = (sl*L.v + sh*R.v)/(sl+sh);
    const double Hbl = (sl*Hl + sh*Hh)/(sl+sh);
    double abl = std::sqrt((cc::gamma-1) * (Hbl - 0.5 * (ubl*ubl + vbl*vbl)));
    double length = std::sqrt(face.nor.x*face.nor.x + face.nor.y*face.nor.y);
    double nx = face.nor.x/length;double ny = face.nor.y/length;
    double unbl = ubl * nx + vbl * ny;
    double utbl = vbl * nx - ubl * ny;

    // 特征值
    double lambda[4] = {std::abs(unbl - abl),std::abs(unbl),std::abs(unbl),std::abs(unbl + abl)}; 

    // 波强度
    const double dp = R.p - L.p;
    const double drho = R.rho - L.rho;
    const double dun = nx * (R.u - L.u) + ny * (R.v - L.v);
    const double dut = nx * (R.v - L.v) - ny * (R.u - L.u);
    double alpha[4] = {(dp-rhobl*abl*dun)/(2*abl*abl),drho-dp/(abl*abl),rhobl*dut,(dp+rhobl*abl*dun)/(2*abl*abl)};
    
    // 特征向量
    double tz[4][4] = {{1,ubl-abl*nx,vbl-abl*ny,Hbl-unbl*abl},
                       {1,ubl,       vbl,       (ubl*ubl+vbl*vbl)*0.5},
                       {0,-ny,       nx,        utbl},
                       {1,ubl+abl*nx,vbl+abl*ny,Hbl+unbl*abl}};

    // 左右无粘通量
    const double unl = L.u * nx + L.v * ny;
    const double unh = R.u * nx + R.v * ny;
    const double FL[4] = {L.rho*unl,L.rho*L.u*unl+L.p*nx,
                          L.rho*L.v*unl+L.p*ny,unl*(L.rho*L.e+L.p)};
    const double FR[4] = {R.rho*unh,R.rho*R.u*unh+R.p*nx,
                          R.rho*R.v*unh+R.p*ny,unh*(R.rho*R.e+R.p)};

    // 形成面上对流通量(最后乘面长)
    for(int i=0;i<4;i++){
        face.convect[i] = 0.5*(FL[i]+FR[i]);
        for(int j=0;j<4;j++){
            face.convect[i] -= 0.5*lambda[j]*alpha[j]*tz[j][i];
        }
        face.convect[i] *= length;
    }
}