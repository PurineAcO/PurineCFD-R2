#pragma once

#include "classconfig.hpp"
#include "config.hpp"
#include <cstdio>
#include <cstdlib>

/*
HALO 虚网格一般只能适用于结构化网格求解,为了满足精度此处默认取3个.
HALO 网格独立存储在cc::GhostList中,从横行0~3是壁面虚拟网格,3~5是远场虚拟网格.
在已经确定是结构化网格的基础上,建议弃用cc::gotocell,而是使用更加广义的gotoHALO
*/

// 检查网格是否满足结构化条件
bool check_if_structed();
// 建立HALO虚网格
void HALO_structer_mesh();
// 更新虚网格物理量
void update_ghost_field();
// 找到HALO框架下任何网格
cc::cell_class* gotoHALO(int n, int s);

inline bool check_if_structed() {
    for (cc::cell_class& cell : cc::CellList)
        if (cell.index != (cell.n - 1) * structer::S_MAX + cell.s) return false;
    return true;
}

inline cc::cell_class* gotoHALO(int n, int s) {
    if (s > structer::S_MAX || n < -3 || n > structer::N_MAX + 3)
        return nullptr;
    else if (n <= 0)
        return &cc::GhostList[std::abs(n) * structer::S_MAX + s - 1];
    else if (n > structer::N_MAX)
        return &cc::GhostList[(n - structer::N_MAX + 2) * structer::S_MAX + s - 1];
    else
        return &cc::gotocell((n - 1) * structer::S_MAX + s);
}

inline void HALO_structer_mesh() {
    if (!structer::ifstructer) { return; }

    // 给虚网格分配内存
    cc::GhostList.reserve(6 * structer::S_MAX);
    for (int layer = 0; layer < 6; layer++) {
        for (int s = 1; s <= structer::S_MAX; s++) {
            cc::cell_class ghost;
            ghost.index = cc::cell_num + layer * structer::S_MAX + s;
            ghost.s = s;
            if (layer < 3)
                ghost.n = -layer;
            else
                ghost.n = structer::N_MAX + layer - 2;
            ghost.ecnt = 0;
            ghost.east = 0;
            ghost.west = 1;
            ghost.north = 2;
            ghost.south = 3;
            cc::GhostList.push_back(ghost);
        }
    }

    // 链接虚网格
    for (int n = -2; n <= structer::N_MAX + 3; n++) {
        if (n >= 1 && n <= structer::N_MAX) continue;
        for (int s = 1; s <= structer::S_MAX; s++) {
            cc::cell_class* ghost = gotoHALO(n, s);
            ghost->nei[ghost->east] = gotoHALO(n, s == structer::S_MAX ? 1 : s + 1);
            ghost->nei[ghost->west] = gotoHALO(n, s == 1 ? structer::S_MAX : s - 1);
            ghost->nei[ghost->north] = gotoHALO(n + 1, s);
            ghost->nei[ghost->south] = gotoHALO(n - 1, s);
        }
    }

    for (cc::face_class* face : cc::WallFaces) {
        bool i = face->nei[0] ? 1 : 0;
        face->nei[i] = gotoHALO(0, face->nei[!i]->s);
        face->nei[!i]->nei[face->nei[!i]->south] = gotoHALO(0, face->nei[!i]->s);
    }
    for (cc::face_class* face : cc::FarFaces) {
        bool i = face->nei[0] ? 1 : 0;
        face->nei[i] = gotoHALO(structer::N_MAX + 1, face->nei[!i]->s);
        face->nei[!i]->nei[face->nei[!i]->north] = gotoHALO(structer::N_MAX + 1, face->nei[!i]->s);
    }

    update_ghost_field();
    printf("HALO: 6 layers, %d ghost cells\n", 6 * structer::S_MAX);
}

inline void update_ghost_field() {

    // 壁面
    for (int n = 0; n >= -2; n--) {
        for (int s = 1; s <= structer::S_MAX; s++) {
            cc::cell_class& ghost = *gotoHALO(n, s);
            const cc::cell_class& inner = *gotoHALO(1 - n, s);
            ghost.phy = inner.phy;
            ghost.otphy = inner.otphy;
            ghost.phy.u = -inner.phy.u;
            ghost.phy.v = -inner.phy.v;
            ghost.tur.miubl = inner.tur.miubl;
        }
    }

    // 压力远场
    for (int n = structer::N_MAX + 1; n <= structer::N_MAX + 3; n++) {
        for (int s = 1; s <= structer::S_MAX; s++) {
            cc::cell_class& ghost = *gotoHALO(n, s);
            const cc::cell_class& inner = *gotoHALO(structer::N_MAX, s);
            cc::face_class* reface = inner.northf;
            ghost.phy = reface->phy;
            ghost.otphy = reface->otphy;
            ghost.tur.miubl = reface->tur.miubl;
        }
    }

    // 耗散项会读邻居的守恒量
    for (cc::cell_class& ghost : cc::GhostList)
        ghost.form_conservative();
}
