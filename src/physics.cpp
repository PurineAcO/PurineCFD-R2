#include "physics.h"
#include <cmath>

double sound_speed(double T) {
  return sqrt(cfd::gamma * cfd::R * T);
}

double sutherland::dynamic_viscosity(double T) {
  return sutherland::mu0 * (T / sutherland::T0) * std::sqrt(T / sutherland::T0) *
         (sutherland::T0 + sutherland::Ts) / (T + sutherland::Ts);
}
