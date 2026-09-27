#pragma once

#include "classconfig.hpp"
#include "config.hpp"

/*
MUSCL只能建立邻接面的左物理量和右物理量,无论是结构网格还是非结构网格都显式依赖了梯度
面上梯度的处理使用了正交和非正交分裂后修正的技术
*/

// 面上中心差分插值
void interpolate_mid(cc::face_class* face);
// MUSCL左右插值
void muscl(cc::face_class* face);
// 面上梯度的处理
void grad_onface(cc::face_class *face);

inline void interpolate_mid(cc::face_class *face){face->face_physic_mid();}

inline void muscl(cc::face_class *face){
    if(face->type != cc::INTER){return;}
    face->phynei[0] = face->nei[0]->phy + face->nei[0]->phgrad*(face->mid - face->nei[0]->center);
    face->phynei[1] = face->nei[1]->phy + face->nei[1]->phgrad*(face->mid - face->nei[1]->center);
}

inline void grad_onface(cc::face_class *face){
    face->phgrad.clear();
    if(face->type != cc::INTER) {
        cc::cell_class* inner = (face->nei[0]->index<=cc::cell_num) ? face->nei[0] : face->nei[1];
        face->phgrad = (face->phy - inner->phy)*(1.0/(face->mid-inner->center).norm());
    }
    cc::vec2 to = face->nei[1]->center - face->nei[0]->center;
    cc::vec2 nor = face->nor * face->outer;
    // 正交部分
    face->phgrad += (face->nei[1]->phy-face->nei[0]->phy) * cc::dot(nor, to*(1.0/to.norm()));
    // 非正交修正
    cc::vec2 weight ={cc::dot(face->nei[1]->center-face->mid,nor),cc::dot(face->mid-face->nei[0]->center,nor)};
    weight *= 1/(weight.x+weight.y);
    cc::vecgrad no_grad = weight.x * face->nei[0]->phgrad + weight.y * face->nei[1]->phgrad;
    cc::vec2 after = nor - cc::dot(nor, to*(1.0/to.norm())) * to*(1.0/to.norm());
    face->phgrad += no_grad * after;
}