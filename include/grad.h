#pragma once
#include "classconfig.h"

// GGCB梯度
void green_gauss_cell_based(cc::cell_class& cell);
// 面上梯度: 内部面取平均, 壁面修正法向分量
void face_gradient(cc::face_class& face);
