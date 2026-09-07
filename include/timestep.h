#pragma once
#include "mesh.h"
#include <array>

namespace rk {
inline constexpr std::array<double, 5> coefficients = {0.25, 1.0 / 6, 0.375, 0.5, 1.0};
}
void compute_local_timestep(cfd::Cell& cell);
