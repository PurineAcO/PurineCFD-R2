#pragma once

#include "classconfig.hpp"
#include "config.hpp"

/*
MUSCL只能建立邻接面的左物理量和右物理量,无论是结构网格还是非结构网格都显式依赖了梯度
面上梯度按照OpenFOAM的fvc::interpolate::dotInterpolate进行重构,边界面上梯度暂使用一阶外推.
*/

// 面上中心差分插值
void interpolate_mid(cc::face_class* face);
// MUSCL左右插值
void muscl(cc::face_class* face);
// 面上梯度的处理
void grad_onface(cc::cell_class &cell);

inline void interpolate_mid(cc::face_class *face){face->face_physic_mid();face->form_otherphy();}

// Barth-Jespersen 限制器: 重构值不超出邻居极值, 保证正性
inline cc::vecp limited_phy(const cc::cell_class* cell,const cc::vec2& off){
    const double c[4] = {cell->phy.rho,cell->phy.u,cell->phy.v,cell->phy.T};
    double lo[4] = {c[0],c[1],c[2],c[3]},hi[4] = {c[0],c[1],c[2],c[3]};
    for(int i=0;i<cell->ecnt;i++){
        const cc::face_class* f = cell->faces[i];
        const cc::cell_class* o = (f->nei[0] == cell) ? f->nei[1] : f->nei[0];
        if(o == nullptr){ continue; }
        const double n[4] = {o->phy.rho,o->phy.u,o->phy.v,o->phy.T};
        for(int k=0;k<4;k++){
            if(n[k] < lo[k]){ lo[k] = n[k]; }
            if(n[k] > hi[k]){ hi[k] = n[k]; }
        }
    }
    const cc::vecp d = cell->phgrad*off;
    const double dk[4] = {d.rho,d.u,d.v,d.T};
    double r[4] = {c[0],c[1],c[2],c[3]};
    for(int k=0;k<4;k++){
        double psi = 1.0;
        if(dk[k] > 0.0){
            const double m = (hi[k] - c[k])/dk[k];
            psi = m < 1.0 ? m : 1.0;
        }else if(dk[k] < 0.0){
            const double m = (lo[k] - c[k])/dk[k];
            psi = m < 1.0 ? m : 1.0;
        }
        psi = psi > 0.0 ? psi : 0.0;
        r[k] = c[k] + psi*dk[k];
    }
    return cc::vecp(r[0],r[1],r[2],r[3]);
}

inline void muscl(cc::face_class *face){
    if(face->type == cc::INTER){
        face->phynei[0] = limited_phy(face->nei[0],face->mid - face->nei[0]->center);
        face->phynei[1] = limited_phy(face->nei[1],face->mid - face->nei[1]->center);
    }else{
        // 边界面只有一侧有格子, 取边界条件给出的面值
        face->phynei[0] = face->phy;
        face->phynei[1] = face->phy;
    }
    // convect_ROE 要用左右引申物理量
    face->otnei[0].form_otphy(face->phynei[0]);
    face->otnei[1].form_otphy(face->phynei[1]);
}

inline void grad_onface(cc::cell_class &cell){
    allface(cell){
        cc::face_class* face = cell.faces[i];
        // 先算局部值再整体赋值: 并行时两侧格子会同时写同一个面
        cc::vecgrad pg;cc::vec2 mg;
        pg.clear();mg.clear();
        if(face->type == cc::INTER){
            // 面梯度只由两侧格子决定, 与调用它的格子无关
            const cc::cell_class* lo = face->nei[0];
            const cc::cell_class* hi = face->nei[1];
            double weight = cc::dot(face->nor,face->mid-lo->center)/
            (cc::dot(face->nor,hi->center-face->mid) + cc::dot(face->nor,face->mid-lo->center));
            pg = weight * lo->phgrad + (1-weight) * hi->phgrad;
            mg = weight * lo->tur.miublgrad + (1-weight) * hi->tur.miublgrad;
        }
        else if(face->type == cc::WALL){
            // 法向导数用一阶壁面律: 格子中心的最小二乘梯度看不到壁面那一跳
            const double y1 = cell.tur.sad;
            const double dx = face->mid.x - cell.center.x;
            const double dy = face->mid.y - cell.center.y;
            const double len = std::sqrt(dx*dx + dy*dy);
            const double nx = dx/len,ny = dy/len,tx = -ny,ty = nx;
            const double du_dn = (face->phy.u - cell.phy.u)/y1;
            const double dv_dn = (face->phy.v - cell.phy.v)/y1;
            const double dT_dn = (face->phy.T - cell.phy.T)/y1;
            const double db_dn = (face->tur.miubl - cell.tur.miubl)/y1;
            const double du_dt = cell.phgrad.ugrad.x*tx + cell.phgrad.ugrad.y*ty;
            const double dv_dt = cell.phgrad.vgrad.x*tx + cell.phgrad.vgrad.y*ty;
            const double dT_dt = cell.phgrad.Tgrad.x*tx + cell.phgrad.Tgrad.y*ty;
            const double db_dt = cell.tur.miublgrad.x*tx + cell.tur.miublgrad.y*ty;
            pg.ugrad = {du_dn*nx + du_dt*tx,du_dn*ny + du_dt*ty};
            pg.vgrad = {dv_dn*nx + dv_dt*tx,dv_dn*ny + dv_dt*ty};
            pg.Tgrad = {dT_dn*nx + dT_dt*tx,dT_dn*ny + dT_dt*ty};
            mg = {db_dn*nx + db_dt*tx,db_dn*ny + db_dt*ty};
        }
        else{
            pg = cell.phgrad;mg = cell.tur.miublgrad;
        }
        face->phgrad = pg;face->tur.miublgrad = mg;
    }
}