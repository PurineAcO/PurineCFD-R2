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
void grad_onface(cc::face_class *face);

inline void interpolate_mid(cc::face_class *face){face->face_physic_mid();}

inline void muscl_face(cc::face_class *face){
    if(face->type != cc::INTER){return;}
    face->phynei[0] = face->nei[0]->phy + face->nei[0]->phgrad*(face->mid - face->nei[0]->center);
    face->phynei[1] = face->nei[1]->phy + face->nei[1]->phgrad*(face->mid - face->nei[1]->center);
}

inline void grad_onface(cc::cell_class *cell){
    for(int i=1;i<cell->ecnt;i++){
        cell->faces[i]->phgrad.clear();cell->faces[i]->tur.miublgrad.clear();
        if(cell->faces[i]->type == cc::INTER){
        // cc::vec2 df = cell->nei[i]->center - cell->center;
        // cc::vec2 nf = cell->faces[i]->nor * ((2*cell->fnorm[i] - 1) *1.0/(cell->faces[i]->nor.norm()));
        // cc::vec2 kf = nf - 1.0/(cc::dot(nf,df)) * df;
        double weight = cc::dot(cell->faces[i]->nor,cell->faces[i]->mid-cell->center)/
        (cc::dot(cell->faces[i]->nor,cell->nei[i]->center-cell->faces[i]->mid) + cc::dot(cell->faces[i]->nor,cell->faces[i]->mid-cell->center));
        cell->faces[i]->phgrad = weight * cell->phgrad + (1-weight) * cell->nei[i]->phgrad;
        cell->faces[i]->tur.miublgrad = weight * cell->tur.miublgrad + (1-weight) * cell->nei[i]->tur.miublgrad;
        }
        else{
        cell->faces[i]->phgrad = cell->phgrad;cell->faces[i]->tur.miublgrad = cell->tur.miublgrad;
        if(cell->faces[i]->type == cc::WALL){cell->faces[i]->phgrad.rhograd.clear();cell->faces[i]->phgrad.Tgrad.clear();}
        }
    }
}