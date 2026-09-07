#include "spalart_allmaras.h"
#include "config.h"
#include "mesh.h"
#include "physics.h"
#include <cmath>

namespace sa {
namespace {

// χ = max(ν̃, 0)/ν，其中 ν = μ/ρ。
double viscosity_ratio(double rho, double mu, double nu_tilde) {
  return rho * (nu_tilde > 0.0 ? nu_tilde : 0.0) / mu;
}
// fv1 将工作变量转换为湍流运动黏度：νt = ν̃*fv1。
double compute_fv1(double chi) {
  return (chi * chi * chi) / (chi * chi * chi + sa::Cv1 * sa::Cv1 * sa::Cv1);
}
// ft2 修正 SA 源项中的产生项和破坏项。
double compute_ft2(double chi) {
  return sa::Ct3 * std::exp(-sa::Ct4 * chi * chi);
}
// fv2 用于修正源项中的涡量尺度。
double compute_fv2(double chi) {
  return 1 - chi / (1 + chi * compute_fv1(chi));
}
// g 是壁面破坏函数 fw 的中间变量。
double compute_g(double r) {
  return r + sa::Cw2 * (r * r * r * r * r * r - r);
}
// fw 控制壁面破坏项随无量纲距离参数 r 的变化。
double compute_fw(double g) {
  constexpr double cw3_squared = sa::Cw3 * sa::Cw3;
  constexpr double cw3_sixth = cw3_squared * cw3_squared * cw3_squared;
  const double g_squared = g * g;
  const double g_sixth = g_squared * g_squared * g_squared;
  return g * std::pow((1 + cw3_sixth) / (cw3_sixth + g_sixth), 1.0 / 6);
}

// SA 的体积源项：产生项 - 壁面破坏项 + 梯度平方项（此处含密度因子）。
double source_term(const cfd::Cell& cell) {
  double mu = sutherland::dynamic_viscosity(cell.flow.T);
  double chi = viscosity_ratio(cell.flow.rho, mu, cell.turbulence.nu_tilde);
  double ft2 = compute_ft2(chi);
  double fv2 = compute_fv2(chi);
  double vorticity = std::abs(cell.flow.vgrad.x - cell.flow.ugrad.y);
  const double scaled_nu =
      cell.turbulence.nu_tilde * cell.turbulence.inverse_wall_distance_squared / (kappa * kappa);
  double modified_vorticity =
      std::max(vorticity + fv2 * scaled_nu, std::max(0.3 * vorticity, 1e-20));
  double production =
      Cb1 * (1 - ft2) * modified_vorticity * cell.flow.rho * cell.turbulence.nu_tilde;
  double r = std::min(scaled_nu / modified_vorticity, rmax);
  double g = compute_g(r);
  double destruction = cell.flow.rho * (Cw1 * compute_fw(g) - Cb1 / kappa / kappa * ft2) *
                       cell.turbulence.nu_tilde * cell.turbulence.nu_tilde *
                       cell.turbulence.inverse_wall_distance_squared;
  double gradient_source =
      Cb2 * inv_sigma * cell.flow.rho *
      cfd::dot(cell.turbulence.nu_tilde_gradient, cell.turbulence.nu_tilde_gradient);
  return production - destruction + gradient_source;
}

} // namespace

double freestream_nu_tilde() {
  if (!config::uses_sa())
    return 0.0;
  const double rho = cfd::freestream.p / (cfd::R * cfd::freestream.T);
  return 3.0 * sutherland::dynamic_viscosity(cfd::freestream.T) / rho;
}

double eddy_viscosity(const cfd::Face& face, double mu) {
  const double chi = viscosity_ratio(face.flow.rho, mu, face.turbulence.nu_tilde);
  return face.flow.rho * compute_fv1(chi) * face.turbulence.nu_tilde;
}

void prepare_face_flux(cfd::Face& face) {
  const double mu = sutherland::dynamic_viscosity(face.flow.T);
  face.sa_diffusivity = face.turbulence.nu_tilde + mu / face.flow.rho;
  face.sa_gradient_flux = face.turbulence.nu_tilde_gradient.x * face.area_normal.x +
                          face.turbulence.nu_tilde_gradient.y * face.area_normal.y;
}

void advance_turbulence(cfd::Cell& cell, double coefficient) {
  double rhs = 0.0;
  double velocity_divergence = 0.0;
  for (int j = 0; j < cell.face_count; j++) {
    cfd::Face* face = cell.faces[j];
    int outward_sign = 2 * cell.normal_points_outward[j] - 1;
    // 一阶迎风取值：流出时取本单元 ν̃，流入时取相邻单元或边界面的 ν̃。
    double outward_volume_flux = outward_sign * face->volume_flux;
    double upwind_nu = cell.turbulence.nu_tilde;
    if (outward_volume_flux < 0.0) {
      if (const cfd::Cell* other = cell.neighbors[j]) {
        upwind_nu = other->turbulence.nu_tilde;
      } else {
        upwind_nu = face->turbulence.nu_tilde;
      }
    }
    double diffusivity = face->sa_diffusivity;
    double outward_gradient_flux = outward_sign * face->sa_gradient_flux;
    rhs += (outward_volume_flux * upwind_nu - inv_sigma * diffusivity * outward_gradient_flux) *
           cell.inverse_volume;
    velocity_divergence += outward_volume_flux * cell.inverse_volume;
  }
  // 非守恒形式的 ν̃ 方程需要扣除 ν̃*div(u)；体积源项除以密度。
  rhs -= source_term(cell) / cell.flow.rho + cell.turbulence.nu_tilde * velocity_divergence;
  double next_nu =
      cell.turbulence.nu_tilde_previous - sa::relaxation_factor * coefficient * cell.local_dt * rhs;
  // 保存下一阶段值，待所有单元右端项算完后由 solver.cpp 统一写回。
  cell.turbulence.nu_tilde_next = std::isfinite(next_nu) ? std::max(next_nu, 0.0) : next_nu;
}

} // namespace sa
