#include "initialization.h"
#include "mesh.h"
#include "physics.h"
#include "spalart_allmaras.h"
#include <cstdio>

void initialize_freestream() {
  const auto& fs = cfd::freestream;
  const double nu_tilde_inf = sa::freestream_nu_tilde();
  for (cfd::Cell& cell : cfd::cells) {
    cell.flow.T = fs.T;
    cell.flow.p = fs.p;
    cell.flow.u = fs.u;
    cell.flow.v = fs.v;
    cell.flow.rho = cell.flow.p / cell.flow.T / cfd::R;
    cell.flow.a = sound_speed(cell.flow.T);
    cell.flow.e =
        cfd::Cv * cell.flow.T + 0.5 * (cell.flow.u * cell.flow.u + cell.flow.v * cell.flow.v);
    cell.turbulence.nu_tilde = nu_tilde_inf;
  }
  for (cfd::Face& face : cfd::faces) {
    face.turbulence.nu_tilde = (face.type == cfd::BoundaryType::wall) ? 0.0 : nu_tilde_inf;
    if (face.type == cfd::BoundaryType::interior) {
      face.interpolate_primitives();
    }
  }
  printf("STD Initialization OK!, u is: %f\n", cfd::cells[0].flow.u);
}
