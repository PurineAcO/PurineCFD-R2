#pragma once

#include "classconfig.hpp"
#include "config.hpp"
#include "dissipation.hpp"
#include <cmath>

/*
Roe 通量对方向是有要求的,也就是说必须保证面上的法向量是L→R,本代码中要求L→R是face.nei[0]→[1]
也就是说在interpolate中已经要求0和1必须完全匹配.face.outer记载了face.nor和face.nei[0]→[1]的方向关系
*/

// 无粘通量
void convect_JST(cc::face_class& face);
// 汇总单元各面的无粘和黏性通量
void assemble_flux(cc::cell_class& cell, char fluxtype = 'R');
// Roe无粘通量
void convect_ROE(cc::face_class& face);

inline void convect_JST(cc::face_class& face) {
    const double rho = face.phy.rho, u = face.phy.u, v = face.phy.v;
    const double H = cc::Cp * face.phy.T + 0.5 * (u * u + v * v);
    const double p = cc::R * rho * face.phy.T;
    const cc::vec4 F = {rho * u, rho * u * u + p, rho * u * v, rho * u * H};
    const cc::vec4 G = {rho * v, rho * u * v, rho * v * v + p, rho * v * H};
    face.convect = face.toface_jacobi(F, G);
    jst::dissipation(face);
}

inline void assemble_flux(cc::cell_class& cell, char fluxtype) {
    cell.convect.clear();
    for (int i = 0; i < cell.ecnt; i++) {
        // TODO:需要在config.json保留一个JST或者ROE的选项,这里准备硬编码
        if (fluxtype == 'J') { cell.convect += (2 * cell.fnorm[i] - 1) * cell.faces[i]->convect; }
        if (fluxtype == 'R') {
            bool islow = (cell.faces[i]->nei[0] == &cell) ? true : false;
            cell.convect += (2 * islow - 1) * cell.faces[i]->convect;
        }
    }
}

inline void convect_ROE(cc::face_class& face) {
    // MUSCL重构后的左右状态
    const cc::vecp& L = face.phynei[0];
    const cc::vecp& R = face.phynei[1];
    const cc::otphy& Lo = face.otnei[0];
    const cc::otphy& Ro = face.otnei[1];
    // Roe平均值
    const double sl = std::sqrt(L.rho), sh = std::sqrt(R.rho);
    const double Hl = cc::Cp * L.T + 0.5 * (L.u * L.u + L.v * L.v);
    const double Hh = cc::Cp * R.T + 0.5 * (R.u * R.u + R.v * R.v);
    const double rhobl = sl * sh;
    const double ubl = (sl * L.u + sh * R.u) / (sl + sh);
    const double vbl = (sl * L.v + sh * R.v) / (sl + sh);
    const double Hbl = (sl * Hl + sh * Hh) / (sl + sh);
    double abl = std::sqrt((cc::gamma - 1) * (Hbl - 0.5 * (ubl * ubl + vbl * vbl)));
    double length = std::sqrt(face.nor.x * face.nor.x + face.nor.y * face.nor.y);
    // 边界面 nei[0] 可能为空, 此时 outer 仍是默认值, 必须现场判向
    const double sign = face.nei[0] != nullptr
                            ? (cc::dot(face.nor, face.mid - face.nei[0]->center) > 0 ? 1.0 : -1.0)
                            : (cc::dot(face.nor, face.nei[1]->center - face.mid) > 0 ? 1.0 : -1.0);
    double nx = face.nor.x / length * sign;
    double ny = face.nor.y / length * sign;
    double unbl = ubl * nx + vbl * ny;
    double utbl = vbl * nx - ubl * ny;

    // 特征值
    double lambda[4] = {std::abs(unbl - abl), std::abs(unbl), std::abs(unbl), std::abs(unbl + abl)};

    // 波强度
    const double dp = Ro.p - Lo.p;
    const double drho = R.rho - L.rho;
    const double dun = nx * (R.u - L.u) + ny * (R.v - L.v);
    const double dut = nx * (R.v - L.v) - ny * (R.u - L.u);
    double alpha[4] = {(dp - rhobl * abl * dun) / (2 * abl * abl), drho - dp / (abl * abl),
                       rhobl * dut, (dp + rhobl * abl * dun) / (2 * abl * abl)};

    // 特征向量
    cc::vec4 tz[4] = {{1, ubl - abl * nx, vbl - abl * ny, Hbl - unbl * abl},
                      {1, ubl, vbl, (ubl * ubl + vbl * vbl) * 0.5},
                      {0, -ny, nx, utbl},
                      {1, ubl + abl * nx, vbl + abl * ny, Hbl + unbl * abl}};

    // 左右无粘通量
    const double unl = L.u * nx + L.v * ny;
    const double unh = R.u * nx + R.v * ny;
    const cc::vec4 FL = {L.rho * unl, L.rho * L.u * unl + Lo.p * nx, L.rho * L.v * unl + Lo.p * ny,
                         unl * (L.rho * Lo.e + Lo.p)};
    const cc::vec4 FR = {R.rho * unh, R.rho * R.u * unh + Ro.p * nx, R.rho * R.v * unh + Ro.p * ny,
                         unh * (R.rho * Ro.e + Ro.p)};

    // 形成面上对流通量(最后乘面长)
    face.convect = 0.5 * (FL + FR);
    for (int i = 0; i < 4; i++) {
        face.convect -= 0.5 * lambda[i] * alpha[i] * tz[i];
    }
    face.convect = face.convect * length;
}
