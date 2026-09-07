#include "timestep.h"
#include "config.h"
#include "mesh.h"
#include "physics.h"
#include "spalart_allmaras.h"
#include <cmath>

// Δt = CFL /（对流谱半径 + 扩散稳定性界 + SA 源项稳定性界）。
void compute_local_timestep(cfd::Cell& cell) {
  const double dsx = cell.projected_face_sum.x;
  const double dsy = cell.projected_face_sum.y;
  double convective_x = (std::abs(cell.flow.u) + cell.flow.a) * dsx;
  double convective_y = (std::abs(cell.flow.v) + cell.flow.a) * dsy;
  double viscous_bound = 0.0;
  {
    double diffusivity_bound =
        sutherland::dynamic_viscosity(cell.flow.T) / cell.flow.rho + cell.turbulence.nu_tilde;
    viscous_bound =
        4.0 * diffusivity_bound * (dsx * dsx + dsy * dsy) * cell.inverse_volume / cfd::Pr;
  }
  double sa_source_bound =
      2.0 * sa::Cw1 * cell.turbulence.nu_tilde * cell.turbulence.inverse_wall_distance_squared;
  cell.local_dt =
      config::settings.cfl /
      ((convective_x + convective_y + viscous_bound) * cell.inverse_volume + sa_source_bound);
}
