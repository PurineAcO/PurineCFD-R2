#include "solver.h"
#include "boundary.h"
#include "config.h"
#include "convergence.h"
#include "fluxes.h"
#include "geometry.h"
#include "gradients.h"
#include "initialization.h"
#include "io.h"
#include "jst.h"
#include "mesh.h"
#include "mesh_reader.h"
#include "spalart_allmaras.h"
#include "timestep.h"
#include "viscous.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <omp.h>
#include <stdexcept>
#include <string>

namespace {
template <typename Action> void for_each_cell(Action action) {
#pragma omp parallel for schedule(static)
  for (auto& cell : cfd::cells)
    action(cell);
}

template <typename Action> void for_each_face(Action action) {
#pragma omp parallel for schedule(static)
  for (auto& face : cfd::faces)
    action(face);
}

void recover_flow() {
  for_each_cell([](auto& cell) {
    cell.recover_primitives();
    cell.update_thermodynamics();
  });
}

void check_flow(int step) {
  int first_invalid = cfd::cell_count;
#pragma omp parallel for schedule(static) reduction(min : first_invalid)
  for (int i = 0; i < cfd::cell_count; ++i) {
    const auto& cell = cfd::cells[i];
    const auto& state = cell.flow;
    if (!std::isfinite(state.rho) || state.rho <= 0.0 || !std::isfinite(state.T) ||
        state.T <= 0.0 || !std::isfinite(state.u) || !std::isfinite(state.v) ||
        !std::isfinite(state.p) || !std::isfinite(state.e) ||
        !std::isfinite(cell.turbulence.nu_tilde) || cell.turbulence.nu_tilde < 0.0) {
      first_invalid = std::min(first_invalid, i);
    }
  }
  if (first_invalid < cfd::cell_count)
    throw std::runtime_error("Invalid flow state at step " + std::to_string(step) + ", cell " +
                             std::to_string(cfd::cells[first_invalid].index));
}

// Q^(k) = Qⁿ - αk*Δt*R(Q^(k-1))；每一阶段都以本步起点 Qⁿ 为基准。
void advance_stage(double coefficient) {
  // 1. 从守恒量恢复单元状态，施加边界条件，再插值内部面的原始变量。
  recover_flow();
  apply_wall_boundary();
  apply_farfield_boundary();
  for_each_face([](auto& face) {
    face.interpolate_primitives();
    face.update_thermodynamics();
    if (face.type == cfd::BoundaryType::interior) {
      face.turbulence.nu_tilde =
          0.5 * (face.cells[0]->turbulence.nu_tilde + face.cells[1]->turbulence.nu_tilde);
    }
  });
  // 2. 由面值计算单元梯度，并准备 JST 压力传感器和离散 Laplace 量。
  for_each_cell([](auto& cell) {
    compute_cell_gradients(cell);
    jst::compute_pressure_sensor(cell);
    jst::compute_state_laplacian(cell);
  });
  // 3. 等全部单元梯度就绪后计算面通量；每个面只由一个线程写入。
  for_each_face([](auto& face) {
    compute_face_gradients(face);
    prepare_convective_flux(face);
    prepare_viscous_flux(face);
    if (config::uses_sa())
      sa::prepare_face_flux(face);
  });
  // 4. 单元按自己的外法向汇总通量。此处会读取相邻单元，不能提前更新状态。
  for_each_cell([coefficient](auto& cell) {
    // 直接使用已保存的守恒量，避免 Q→原始量→Q 的重复转换。
    assemble_face_fluxes(cell);
    jst::compute_dissipation(cell);
    if (config::uses_sa())
      sa::advance_turbulence(cell, coefficient);
  });
  // 5. 上一遍循环结束后，统一推进守恒量及 SA 工作变量。
  for_each_cell([coefficient](auto& cell) {
    for (int component = 0; component < 4; ++component) {
      const double residual =
          (cell.convective_flux[component] - cell.dissipation_terms.flux[component] -
           cell.viscous_flux[component]) *
          cell.inverse_volume;
      cell.conservative[component] =
          cell.previous_conservative[component] - coefficient * cell.local_dt * residual;
    }
    if (config::uses_sa())
      cell.turbulence.nu_tilde = cell.turbulence.nu_tilde_next;
  });
}

} // namespace

void solve_steady_flow() {
  const auto& options = config::settings;
  read_mesh(options.mesh_path);
  initialize_geometry();
  if (config::uses_sa())
    for_each_cell(cache_wall_distance);
  initialize_freestream();
  for_each_cell([](auto& cell) { cell.update_conservative(); });
  check_flow(0);
  dump_field(0);

  int step = 0;
  int last_dump = 0;
  double initial_update = -1.0;
  bool converged = false;
  while (step < options.max_steps) {
    ++step;
    for_each_cell([](auto& cell) {
      cell.save_previous_state();
      cell.turbulence.nu_tilde_previous = cell.turbulence.nu_tilde;
      compute_local_timestep(cell);
    });
    // 一个伪时间步包含五个 RK 阶段；local_dt 可因单元而异。
    for (double coefficient : rk::coefficients)
      advance_stage(coefficient);
    recover_flow();
    check_flow(step);
    convergence::report_update(step);
    if (step % options.dump_interval == 0) {
      dump_field(step);
      last_dump = step;
    }
    if (step % options.convergence_interval == 0) {
      const double update = convergence::normalized_max_update();
      if (initial_update < 0.0)
        initial_update = update;
      if (update <= initial_update * 1e-4 && update < 1e-6) {
        std::printf("Converged at step %d, normalized_update=%.6e\n", step, update);
        converged = true;
        break;
      }
    }
  }
  if (last_dump != step)
    dump_field(step);
  if (!converged)
    std::printf("Iteration limit reached; convergence criterion not satisfied.\n");
  std::printf("Total step: %d\n", step);
}
