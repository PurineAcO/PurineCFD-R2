#include "fluxes.h"

void prepare_convective_flux(cfd::Face& face) {
  const double rho = face.flow.rho, u = face.flow.u, v = face.flow.v;
  const double H = cfd::Cp * face.flow.T + 0.5 * (u * u + v * v);
  const double p = cfd::R * rho * face.flow.T;
  const double F[4] = {rho * u, rho * u * u + p, rho * u * v, rho * u * H};
  const double G[4] = {rho * v, rho * u * v, rho * v * v + p, rho * v * H};
  for (int j = 0; j < 4; ++j)
    face.convective_flux[j] = F[j] * face.area_normal.x + G[j] * face.area_normal.y;
}

void assemble_face_fluxes(cfd::Cell& cell) {
  for (int j = 0; j < 4; ++j) {
    cell.convective_flux[j] = 0.0;
    cell.viscous_flux[j] = 0.0;
  }
  // 面法向固定；用本单元的方向符号转换为外向通量，按四个面的固定顺序求和。
  for (int i = 0; i < cell.face_count; ++i) {
    const auto& face = *cell.faces[i];
    const int outer = 2 * cell.normal_points_outward[i] - 1;
    for (int j = 0; j < 4; ++j) {
      cell.convective_flux[j] += outer * face.convective_flux[j];
      cell.viscous_flux[j] += outer * face.viscous_flux[j];
    }
  }
}
