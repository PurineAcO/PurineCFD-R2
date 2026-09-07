#include "jst.h"
#include "mesh.h"
#include <cmath>

void jst::compute_pressure_sensor(cfd::Cell& cell) {
  double up = 0, down = 0;
  for (int i = 0; i < cell.face_count; i++) {
    cfd::Cell* neighbor = cell.neighbors[i];
    if (neighbor == nullptr) {
      continue;
    }
    up += std::abs(neighbor->flow.p - cell.flow.p);
    down += neighbor->flow.p + cell.flow.p;
  }
  cell.dissipation_terms.pressure_sensor = (down > 0) ? up / down : 0.0;
}

void jst::compute_state_laplacian(cfd::Cell& cell) {
  for (int j = 0; j < 4; j++) {
    cell.dissipation_terms.laplacian[j] = 0.0;
  }
  for (int i = 0; i < cell.face_count; i++) {
    cfd::Cell* neighbor = cell.neighbors[i];
    if (neighbor == nullptr) {
      continue;
    }
    for (int j = 0; j < 4; j++) {
      cell.dissipation_terms.laplacian[j] += neighbor->conservative[j] - cell.conservative[j];
    }
  }
}

void jst::compute_dissipation(cfd::Cell& cell) {
  for (int j = 0; j < 4; j++) {
    cell.dissipation_terms.flux[j] = 0.0;
  }
  for (int i = 0; i < cell.face_count; i++) {
    cfd::Face* face = cell.faces[i];
    cfd::Cell* neighbor = cell.neighbors[i];
    if (neighbor == nullptr) {
      continue;
    }

    // 计算该面的流动谱半径
    const double lam = face->spectral_radius;
    if (lam < 1e-30) {
      continue;
    }

    // 形成该面上的eps
    double eps2 = jst::k2 * std::max(cell.dissipation_terms.pressure_sensor,
                                     neighbor->dissipation_terms.pressure_sensor);
    double eps4 = std::max(0.0, jst::k4 - eps2);

    // 形成耗散项
    for (int j = 0; j < 4; j++) {
      cell.dissipation_terms.flux[j] +=
          lam * eps2 * (neighbor->conservative[j] - cell.conservative[j]);
      cell.dissipation_terms.flux[j] +=
          lam * eps4 *
          (cell.dissipation_terms.laplacian[j] - neighbor->dissipation_terms.laplacian[j]);
    }
  }
}
