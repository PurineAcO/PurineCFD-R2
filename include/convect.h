#pragma once
#include "classconfig.h"

// 无粘通量
void convect_JST(cc::face_class& face);
// 汇总单元各面的无粘和黏性通量
void assemble_flux(cc::cell_class& cell);
