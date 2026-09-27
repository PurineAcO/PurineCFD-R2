#pragma once

#include "classconfig.hpp"
#include "config.hpp"

/*
MUSCL只能建立邻接面的左物理量和右物理量,无论是结构网格还是非结构网格都显式依赖了梯度
*/

// 面上中心差分插值
void interpolate_mid(cc::face_class* face);
// MUSCL建立面上左右值
void muscl(cc::face_class* face);


inline void interpolate_mid(cc::face_class *face){face->face_physic_mid();}

inline void muscl(cc::face_class *face){
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