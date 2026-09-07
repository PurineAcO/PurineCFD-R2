#include "viscous.h"
#include "config.h"
#include "physics.h"
#include "spalart_allmaras.h"
#include <cmath>

void prepare_viscous_flux(cfd::Face& face) {
  face.volume_flux = face.flow.u * face.area_normal.x + face.flow.v * face.area_normal.y;
  face.spectral_radius = std::abs(face.volume_flux) + face.flow.a * face.length;
  const double mu = sutherland::dynamic_viscosity(face.flow.T);
  const double mut = config::uses_sa() ? sa::eddy_viscosity(face, mu) : 0.0;
  const double mu_eff = mut + mu;
  const double tau_xx = mu_eff * (4.0 / 3 * face.flow.ugrad.x - 2.0 / 3 * face.flow.vgrad.y);
  const double tau_yy = mu_eff * (4.0 / 3 * face.flow.vgrad.y - 2.0 / 3 * face.flow.ugrad.x);
  const double tau_xy = mu_eff * (face.flow.ugrad.y + face.flow.vgrad.x);
  const double lambda_eff = mu / cfd::Pr + mut / sa::Prt;
  const cfd::Vector2 q = {-lambda_eff * cfd::Cp * face.flow.Tgrad.x,
                          -lambda_eff * cfd::Cp * face.flow.Tgrad.y};
  face.viscous_flux[0] = 0.0;
  face.viscous_flux[1] = tau_xx * face.area_normal.x + tau_xy * face.area_normal.y;
  face.viscous_flux[2] = tau_xy * face.area_normal.x + tau_yy * face.area_normal.y;
  face.viscous_flux[3] = (face.flow.u * tau_xx + face.flow.v * tau_xy - q.x) * face.area_normal.x +
                         (face.flow.u * tau_xy + face.flow.v * tau_yy - q.y) * face.area_normal.y;
}
