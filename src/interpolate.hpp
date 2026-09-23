#pragma once

#include "classconfig.hpp"
#include "config.hpp"

// 面上中心差分插值
void interpolate_mid(cc::face_class* face);
// 面上高阶中心差分,只建立左右值
void mid_2nd_lr(cc::face_class* face);


inline void interpolate_mid(cc::face_class *face){face->face_physic_mid();}

inline void mid_2nd_lr(cc::face_class *face){
    // 使用之前,必须已经建立起来梯度,不建立面上的物理量
    if(face->type != cc::INTER){return;}
    face->lowp.rho = face->low->phy.rho + cc::dot(face->low->phy.rhograd, (face->mid - face->low->center));
    face->highp.rho = face->high->phy.rho + cc::dot(face->high->phy.rhograd,face->mid - face->high->center);
    face->lowp.u = face->low->phy.u + cc::dot(face->low->phy.ugrad,(face->mid - face->low->center));
    face->highp.u = face->high->phy.u + cc::dot(face->high->phy.ugrad, (face->mid - face->high->center));
    face->lowp.v = face->low->phy.v + cc::dot(face->low->phy.vgrad,(face->mid - face->low->center));
    face->highp.v = face->high->phy.v + cc::dot(face->high->phy.vgrad, (face->mid - face->high->center));
    face->lowp.T = face->low->phy.T + cc::dot(face->low->phy.Tgrad,(face->mid - face->low->center));
    face->highp.T = face->high->phy.T + cc::dot(face->high->phy.Tgrad, (face->mid - face->high->center));
}