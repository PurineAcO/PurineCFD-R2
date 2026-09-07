#include "boundary.h"
#include "mesh.h"
#include "physics.h"
#include <cmath>

void apply_wall_boundary() {
  for (cfd::Face* wall : cfd::wall_faces) {
    cfd::Cell* c = cfd::boundary_cell(wall);
    wall->flow.u = 0.0;
    wall->flow.v = 0.0;
    wall->flow.T = c->flow.T;
    wall->flow.rho = c->flow.rho;
    wall->turbulence.nu_tilde = 0.0;
  }
}

void apply_farfield_boundary() {
  double rho_inf = cfd::freestream.p / (cfd::R * cfd::freestream.T);
  double a_inf = sound_speed(cfd::freestream.T);
  double nu_inf = 3.0 * sutherland::dynamic_viscosity(cfd::freestream.T) / rho_inf;
  for (cfd::Face* far : cfd::farfield_faces) {

    cfd::Cell* c = cfd::boundary_cell(far);
    double nx = far->area_normal.x, ny = far->area_normal.y;
    double len = std::sqrt(nx * nx + ny * ny);
    nx /= len;
    ny /= len;
    if (nx * (far->midpoint.x - c->center.x) + ny * (far->midpoint.y - c->center.y) < 0) {
      nx = -nx;
      ny = -ny;
    }

    // Riemann 不变量计算
    double a = sound_speed(c->flow.T);
    double vn = c->flow.u * nx + c->flow.v * ny;
    double vt = -1 * c->flow.u * ny + c->flow.v * nx;
    double vn_inf = cfd::freestream.u * nx + cfd::freestream.v * ny;
    double vt_inf = -cfd::freestream.u * ny + cfd::freestream.v * nx;
    double Rp = vn + 2.0 * a / (cfd::gamma - 1.0);         // 第一不变量,出波
    double Rm = vn_inf - 2.0 * a_inf / (cfd::gamma - 1.0); // 第二不变量,入波
    double vn_star = 0.5 * (Rp + Rm);
    double a_star = 0.25 * (cfd::gamma - 1.0) * (Rp - Rm);
    double s, vt_star;
    if (vn_star >= 0.0) {                                // 出流,取内部
      s = c->flow.p / std::pow(c->flow.rho, cfd::gamma); // 第三不变量,熵的衍生物
      vt_star = vt;                                      // 第四不变量,切向速度
    } else {                                             // 入流,取来流
      s = cfd::freestream.p / std::pow(rho_inf, cfd::gamma);
      vt_star = vt_inf;
    }

    // 还原物理量
    far->flow.u = vn_star * nx - vt_star * ny;
    far->flow.v = vn_star * ny + vt_star * nx;
    far->flow.rho = std::pow(a_star * a_star / (cfd::gamma * s), 1.0 / (cfd::gamma - 1.0));
    far->flow.p = s * std::pow(far->flow.rho, cfd::gamma);
    far->flow.T = far->flow.p / (cfd::R * far->flow.rho);
    far->turbulence.nu_tilde = nu_inf;
  }
}
