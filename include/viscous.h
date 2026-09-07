#pragma once
#include "mesh.h"

// 两种模型共用的黏性应力和热通量；层流模式仅使用分子黏度。
void prepare_viscous_flux(cfd::Face& face);
