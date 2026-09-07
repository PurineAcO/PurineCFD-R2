#include "gradients.h"
#include <cmath>

// Green–Gauss：∇φ = (1/V) Σf φf*n_f*Δs_f。
void compute_cell_gradients(cfd::Cell& cell) {
  cell.flow.ugrad = {0.0, 0.0};
  cell.flow.vgrad = {0.0, 0.0};
  cell.flow.Tgrad = {0.0, 0.0};
  cell.turbulence.nu_tilde_gradient = {0.0, 0.0};
  for (int i = 0; i < cell.face_count; i++) {
    const double s = (2 * cell.normal_points_outward[i] - 1) * cell.inverse_volume;
    cell.flow.ugrad += (cell.faces[i]->flow.u * s) * cell.faces[i]->area_normal;
    cell.flow.vgrad += (cell.faces[i]->flow.v * s) * cell.faces[i]->area_normal;
    cell.flow.Tgrad += (cell.faces[i]->flow.T * s) * cell.faces[i]->area_normal;
    cell.turbulence.nu_tilde_gradient +=
        (cell.faces[i]->turbulence.nu_tilde * s) * cell.faces[i]->area_normal;
  }
}

// 内部面平均两侧梯度；壁面修正法向导数，以满足无滑移和绝热条件。
void compute_face_gradients(cfd::Face& face) {
  if (face.type == cfd::BoundaryType::interior) {
    face.flow.ugrad = 0.5 * (face.cells[0]->flow.ugrad + face.cells[1]->flow.ugrad);
    face.flow.vgrad = 0.5 * (face.cells[0]->flow.vgrad + face.cells[1]->flow.vgrad);
    face.flow.Tgrad = 0.5 * (face.cells[0]->flow.Tgrad + face.cells[1]->flow.Tgrad);
    face.turbulence.nu_tilde_gradient = 0.5 * (face.cells[0]->turbulence.nu_tilde_gradient +
                                               face.cells[1]->turbulence.nu_tilde_gradient);
    return;
  }
  face.flow.ugrad = face.flow.vgrad = face.flow.Tgrad =
      face.turbulence.nu_tilde_gradient = {0.0, 0.0};
  if (face.type != cfd::BoundaryType::wall) {
    return;
  }
  const auto& cell = *cfd::boundary_cell(&face);
  double nx = face.area_normal.x, ny = face.area_normal.y;
  double length = std::hypot(nx, ny);
  nx /= length;
  ny /= length;
  double dx = face.midpoint.x - cell.center.x, dy = face.midpoint.y - cell.center.y;
  if (nx * dx + ny * dy < 0.0) {
    nx = -nx;
    ny = -ny;
  }
  double d = nx * dx + ny * dy;
  auto corrected = [&](cfd::Vector2 gradient, double wall, double interior) {
    double correction = (wall - interior) / d - gradient.x * nx - gradient.y * ny;
    return cfd::Vector2{gradient.x + correction * nx, gradient.y + correction * ny};
  };
  face.flow.ugrad = corrected(cell.flow.ugrad, 0.0, cell.flow.u);
  face.flow.vgrad = corrected(cell.flow.vgrad, 0.0, cell.flow.v);
  face.flow.Tgrad = corrected(cell.flow.Tgrad, cell.flow.T, cell.flow.T);
  face.turbulence.nu_tilde_gradient =
      corrected(cell.turbulence.nu_tilde_gradient, 0.0, cell.turbulence.nu_tilde);
}
