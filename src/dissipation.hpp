#pragma once

#include "classconfig.hpp"
#include <cmath>

/*说明:在一些书籍中,要求对JST的dissipation加入几何权theta_IJ,本求解器未做处理.
由于JST不会是我们的最终选择,因此这只会是一个过渡.
*/

namespace jst {
    inline constexpr double k2 = 0.5;       // 二阶阻尼
    inline constexpr double k4 = 1.0/64;    // 四阶阻尼

    // 激波检测器
    void shockwave_recognize(cc::cell_class& cell);
    // 单元Laplace
    void laplace_dissipation(cc::cell_class& cell);
    // 形成JST耗散项
    void JST_dissipation(cc::cell_class& cell);
}

inline void jst::shockwave_recognize(cc::cell_class& cell){
    double up = 0,down = 0;
    for(int i=0;i<cell.ecnt;i++){
        cc::cell_class* neighbor = cell.nei[i];
        if(neighbor == nullptr){
            continue;
        }
        up += std::abs(neighbor->phy.p - cell.phy.p);
        down += neighbor->phy.p + cell.phy.p;
    }
    cell.diss.Y = (down > 0) ? up/down : 0.0;
}

inline void jst::laplace_dissipation(cc::cell_class& cell){
    for(int j=0;j<4;j++){
        cell.diss.L[j] = 0.0;
    }
    for(int i=0;i<cell.ecnt;i++){
        cc::cell_class* neighbor = cell.nei[i];
        if(neighbor == nullptr){
            continue;
        }
        for(int j=0;j<4;j++){
            cell.diss.L[j] += neighbor->conser[j] - cell.conser[j];
        }
    }
}

inline void jst::JST_dissipation(cc::cell_class& cell){
    for(int j=0;j<4;j++){
        cell.diss.Fd[j] = 0.0;
    }
    for(int i=0;i<cell.ecnt;i++){
        cc::face_class* face = cell.faces[i];
        cc::cell_class* neighbor = cell.nei[i];
        if(neighbor == nullptr){
            continue;
        }
        const double lam = face->lam;
        if(lam < 1e-30){
            continue;
        }
        double eps2 = jst::k2*std::max(cell.diss.Y,neighbor->diss.Y);
        double eps4 = std::max(0.0,jst::k4 - eps2);
        for(int j=0;j<4;j++){
            cell.diss.Fd[j] += lam*eps2*(neighbor->conser[j] - cell.conser[j]);
            cell.diss.Fd[j] += lam*eps4*(cell.diss.L[j] - neighbor->diss.L[j]);
        }
    }
}
