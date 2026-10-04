#pragma once

#include "classconfig.hpp"
#include <cmath>

// JST 人工耗散: 二阶压力开关 + 四阶背景耗散, 作用在中心通量上
namespace jst {
inline constexpr double k2 = 0.5;
inline constexpr double k4 = 1.0 / 64;
// 激波(压力)检测器
void shockwave_recognize(cc::cell_class& cell);
// 单元的伪Laplace
void laplace_dissipation(cc::cell_class& cell);
// 面上耗散, 直接从 face.convect 里减掉
void dissipation(cc::face_class& face);
} // namespace jst

inline void jst::shockwave_recognize(cc::cell_class& cell) {
    double up = 0.0, down = 0.0;
    for (int i = 0; i < cell.ecnt; i++) {
        cc::cell_class* neighbor = cell.nei[i];
        if (neighbor == nullptr) { continue; }
        up += std::abs(neighbor->otphy.p - cell.otphy.p);
        down += neighbor->otphy.p + cell.otphy.p;
    }
    cell.diss.Y = (down > 0.0) ? up / down : 0.0;
}

inline void jst::laplace_dissipation(cc::cell_class& cell) {
    double L[4] = {};
    for (int i = 0; i < cell.ecnt; i++) {
        cc::cell_class* neighbor = cell.nei[i];
        if (neighbor == nullptr) { continue; }
        L[0] += neighbor->conser.c - cell.conser.c;
        L[1] += neighbor->conser.x - cell.conser.x;
        L[2] += neighbor->conser.y - cell.conser.y;
        L[3] += neighbor->conser.e - cell.conser.e;
    }
    for (int j = 0; j < 4; j++) {
        cell.diss.L[j] = L[j];
    }
}

inline void jst::dissipation(cc::face_class& face) {
    if (face.type != cc::INTER) { return; }
    const cc::cell_class* lo = face.nei[0];
    const cc::cell_class* hi = face.nei[1];
    // 谱半径带面面积 (|V_n|+a)*S
    const double lam = (std::abs(face.un) + face.otphy.a) * face.nor.norm();
    const double eps2 = k2 * std::max(lo->diss.Y, hi->diss.Y);
    const double eps4 = std::max(0.0, k4 - eps2);
    const double dU[4] = {hi->conser.c - lo->conser.c, hi->conser.x - lo->conser.x,
                          hi->conser.y - lo->conser.y, hi->conser.e - lo->conser.e};
    // d = eps2*dU - eps4*dLap, dLap 必须取 (R 的Laplace - L 的Laplace)
    const double dL[4] = {hi->diss.L[0] - lo->diss.L[0], hi->diss.L[1] - lo->diss.L[1],
                          hi->diss.L[2] - lo->diss.L[2], hi->diss.L[3] - lo->diss.L[3]};
    // 耗散定义在 nei[0]->nei[1] 方向, 这里连同 nor 一起翻号
    const double k = lam * (2 * face.outer - 1);
    face.convect.c -= k * (eps2 * dU[0] - eps4 * dL[0]);
    face.convect.x -= k * (eps2 * dU[1] - eps4 * dL[1]);
    face.convect.y -= k * (eps2 * dU[2] - eps4 * dL[2]);
    face.convect.e -= k * (eps2 * dU[3] - eps4 * dL[3]);
}
