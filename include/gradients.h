#pragma once
#include "mesh.h"

// 用面值通过 Green–Gauss 公式计算单元梯度。
void compute_cell_gradients(cfd::Cell& cell);
// 插值内部面梯度，并施加壁面法向梯度修正。
void compute_face_gradients(cfd::Face& face);
