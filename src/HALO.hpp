#pragma once

#include "classconfig.hpp"
#include <cstdio>

// 建立HALO虚网格
void HALO_structer_mesh();
// 更新虚网格物理量
void update_ghost_field();
// 找到HALO网格
cc::cell_class& ghost_at(int layer,int s);
// 按外法向定向面的 low/high 侧
void set_face_direction();

// 虚网格不带几何: 环向索引与物理量仅此而已
inline cc::cell_class& ghost_at(int layer,int s){
    return cc::GhostList[layer*structer::S_MAX + s - 1];
}

inline void HALO_structer_mesh(){
    if(!structer::ifstructer){return;}
    const int smax = structer::S_MAX;
    const int nmax = structer::N_MAX;
    const int base = cc::cell_num;
    const int layer_n[6] = {0,-1,-2,nmax+1,nmax+2,nmax+3};

    // 给虚网格分配内存
    cc::GhostList.reserve(6*structer::S_MAX);
    for(int layer=0;layer<6;layer++){
        for(int s=1;s<=smax;s++){
            cc::cell_class ghost;
            ghost.index = base + layer*smax + s;
            ghost.ecnt = 0;
            ghost.s = s;
            ghost.n = layer_n[layer];
            ghost.east = 0;
            ghost.west = 1;
            ghost.north = 2;
            ghost.south = 3;
            cc::GhostList.push_back(ghost);
        }
    }

    // 链接虚网格
    for(int layer=0;layer<6;layer++){
        for(int s=1;s<=smax;s++){
            cc::cell_class& ghost = ghost_at(layer,s);
            ghost.nei[ghost.east] = &ghost_at(layer,s == smax ? 1 : s+1);
            ghost.nei[ghost.west] = &ghost_at(layer,s == 1 ? smax : s-1);
            if(layer == 0){
                ghost.nei[ghost.north] = &cc::CellList[s-1];              // 真实 n = 1
            }else if(layer < 3){
                ghost.nei[ghost.north] = &ghost_at(layer-1,s);
            }else if(layer < 5){
                ghost.nei[ghost.north] = &ghost_at(layer+1,s);
            }
            if(layer > 0 && layer < 3){
                ghost.nei[ghost.south] = &ghost_at(layer+1,s);
            }else if(layer == 3){
                ghost.nei[ghost.south] = &cc::CellList[(nmax-1)*smax + s - 1]; // 真实 n = nmax
            }else if(layer > 3){
                ghost.nei[ghost.south] = &ghost_at(layer-1,s);
            }
        }
    }

    // 边界单元与边界面都接上虚网格
    for(cc::cell_class& cell : cc::CellList){
        for(int i=0;i<cell.ecnt;i++){
            cc::face_class* face = cell.faces[i];
            if(face->type == cc::INTER){
                continue;
            }
            cc::cell_class* ghost = (face->type == cc::WALL) ? &ghost_at(0,cell.s)
                                                             : &ghost_at(3,cell.s);
            cell.nei[i] = ghost;
            if(face->nei[0] == nullptr){
                face->nei[0] = ghost;
            }else{
                face->nei[1] = ghost;
            }
        }
    }

    update_ghost_field();
    printf("HALO: 6 layers, %d ghost cells\n",6*smax);
}

// 面的 low/high 按外法向定向: nor 从 low 指向 high
// 只用格的 fnorm 判定, 不碰任何坐标
inline void set_face_direction(){
    for(cc::cell_class& cell : cc::CellList){
        for(int i=0;i<cell.ecnt;i++){
            cc::cell_class* other = cell.nei[i];
            if(other == nullptr){
                continue;
            }
            cc::face_class* face = cell.faces[i];
            if(cell.fnorm[i]){
                face->low = &cell;
                face->high = other;
            }else{
                face->low = other;
                face->high = &cell;
            }
        }
    }
}

inline void update_ghost_field(){
    const int smax = structer::S_MAX;
    const int nmax = structer::N_MAX;

    // 壁面
    for(int layer=0;layer<3;layer++){
        for(int s=1;s<=smax;s++){
            cc::cell_class& ghost = ghost_at(layer,s);
            const cc::cell_class& inner = cc::CellList[layer*smax + s - 1];
            ghost.phy.rho = inner.phy.rho;
            ghost.phy.p = inner.phy.p;
            ghost.phy.T = inner.phy.T;
            ghost.phy.u = -inner.phy.u;
            ghost.phy.v = -inner.phy.v;
            ghost.phy.a = get_sonic_velocity(ghost.phy.T);
            ghost.phy.e = get_energy(ghost.phy);
            ghost.tur.miubl = -inner.tur.miubl;
        }
    }

    // 压力远场
    const double rho_inf = FAR_DEFINE.p/(cc::R*FAR_DEFINE.T);
    const double miubl_inf = 3.0*sutherland::dynamic_viscosity(FAR_DEFINE.T)/rho_inf;
    for(int layer=3;layer<6;layer++){
        for(int s=1;s<=smax;s++){
            cc::cell_class& ghost = ghost_at(layer,s);
            // 虚网格没有面, 远场面要从最外层真实单元上取
            const cc::cell_class& inner = cc::CellList[(nmax-1)*smax + s - 1];
            cc::face_class* reface = inner.northf;
            ghost.phy.rho = rho_inf;
            ghost.phy.u = reface->phy.u;
            ghost.phy.v = reface->phy.v;
            ghost.phy.T = reface->phy.T;
            ghost.phy.p = reface->phy.p;
            ghost.phy.a = reface->phy.a;
            ghost.phy.e = reface->phy.e;
            ghost.tur.miubl = miubl_inf;
        }
    }

    // 耗散项会读邻居的守恒量
    for(cc::cell_class& ghost : cc::GhostList){
        ghost.form_conservative();
    }
}
