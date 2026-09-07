#include "classconfig.h"
#include "config.h"
#include "grad.h"
#include <cmath>

// 壁面一面: 用壁面格心值替换梯度法向分量为单侧差分,恢复壁面剪切/热流
static constexpr bool kWallShear = false;   // 实验开关: 壁面剪切单侧重构
static cc::vec2 wall_face_grad(cc::vec2 cg, double fval, double cval,
                               double nx, double ny, double invd){
    double cn = cg.x*nx + cg.y*ny;
    double fn = (fval - cval) * invd;
    double dn = fn - cn;
    return {cg.x + dn*nx, cg.y + dn*ny};
}

void interpolate_mid(cc::cell_class &cell){
    for(int i=0;i<cell.ecnt;i++){
        cell.nei[i]->face_physic_mid();
        cell.nei[i]->form_physic();
        if(cell.nei[i]->type == cc::INTER){
            cell.nei[i]->tur.miubl = (cell.nei[i]->nei[0]->tur.miubl + cell.nei[i]->nei[1]->tur.miubl)/2;
        }
    }
    green_gauss_cell_based(cell);
    for(int i=0;i<cell.ecnt;i++){
        cc::face_class* face = cell.nei[i];
        if(face->type == cc::INTER){
            face->phy.ugrad     = 0.5 * (face->nei[0]->phy.ugrad + face->nei[1]->phy.ugrad);
            face->phy.vgrad     = 0.5 * (face->nei[0]->phy.vgrad + face->nei[1]->phy.vgrad);
            face->phy.Tgrad     = 0.5 * (face->nei[0]->phy.Tgrad + face->nei[1]->phy.Tgrad);
            face->tur.miublgrad = 0.5 * (face->nei[0]->tur.miublgrad + face->nei[1]->tur.miublgrad);
        }else if(kWallShear && cc::viscous && face->type == cc::WALL){
            // 无滑移壁面: 沿壁面法向单侧差分重构梯度
            double nx = face->nor.x, ny = face->nor.y;
            double len = std::sqrt(nx*nx + ny*ny);
            double dx = face->mid.x - cell.center.x;
            double dy = face->mid.y - cell.center.y;
            if(nx*dx + ny*dy < 0.0){ nx = -nx; ny = -ny; }
            nx /= len; ny /= len;
            double d = dx*nx + dy*ny;
            if(d <= 0.0){ d = std::sqrt(dx*dx + dy*dy); }
            double invd = 1.0/d;
            face->phy.ugrad     = wall_face_grad(cell.phy.ugrad,     face->phy.u,     cell.phy.u,     nx, ny, invd);
            face->phy.vgrad     = wall_face_grad(cell.phy.vgrad,     face->phy.v,     cell.phy.v,     nx, ny, invd);
            face->phy.Tgrad     = wall_face_grad(cell.phy.Tgrad,     face->phy.T,     cell.phy.T,     nx, ny, invd);
            face->tur.miublgrad = wall_face_grad(cell.tur.miublgrad, face->tur.miubl, cell.tur.miubl, nx, ny, invd);
        }else{
            face->phy.ugrad = {0.0,0.0};
            face->phy.vgrad = {0.0,0.0};
            face->phy.Tgrad = {0.0,0.0};
            face->tur.miublgrad = {0.0,0.0};
        }
    }
}