#include "convergence.h"
#include "mesh.h"
#include "physics.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

namespace convergence {

namespace {
struct Peak {
  double value = -1.0;
  int cell = 0;
};
// 相同更新量取编号较小的单元，使日志不依赖线程调度。
Peak larger(Peak a, Peak b) {
  if (b.value > a.value || (b.value == a.value && b.cell < a.cell))
    return b;
  return a;
}
#pragma omp declare reduction(peak:Peak \
                              : omp_out = larger(omp_out, omp_in)) initializer(omp_priv = {})

constexpr int fields = 5; // rho,u,v,e,nu_tilde；层流的第五分量恒为零。
std::vector<std::array<double, fields>> previous_state;
bool first_report = true;
double max_updates[fields]; // 各原始变量的最大绝对更新量
} // namespace

double normalized_max_update() {
  double rho = cfd::freestream.p / (cfd::R * cfd::freestream.T);
  double U = std::max(std::hypot(cfd::freestream.u, cfd::freestream.v), 1e-20);
  double scales[5] = {rho, U, U, cfd::Cv * cfd::freestream.T,
                      3.0 * sutherland::dynamic_viscosity(cfd::freestream.T) / rho};
  double maximum = 0.0;
  for (int s = 0; s < fields; s++) {
    maximum = std::max(maximum, max_updates[s] / scales[s]);
  }
  return maximum;
}

// 分量最大值用于收敛判断；未归一化的向量模仅用于定位日志中的 maxcell。
void report_update(int step) {
  if (previous_state.empty()) {
    previous_state.resize(cfd::cell_count);
  }
  double max_delta[fields] = {0, 0, 0, 0, 0};
  Peak worst;
#pragma omp parallel for schedule(static) reduction(max \
                                                    : max_delta[:fields]) reduction(peak \
                                                                                    : worst)
  for (int idx = 0; idx < cfd::cell_count; ++idx) {
    const auto& cell = cfd::cells[idx];
    double state[fields] = {cell.flow.rho, cell.flow.u, cell.flow.v, cell.flow.e,
                            cell.turbulence.nu_tilde};
    if (first_report) {
      for (int s = 0; s < fields; s++) {
        previous_state[idx][s] = state[s];
      }
    } else {
      double change[fields];
      double update_norm = 0.0;
      for (int s = 0; s < fields; s++) {
        change[s] = state[s] - previous_state[idx][s];
        update_norm += change[s] * change[s];
      }
      update_norm = std::sqrt(update_norm);
      for (int s = 0; s < fields; s++) {
        if (std::fabs(change[s]) > max_delta[s]) {
          max_delta[s] = std::fabs(change[s]);
        }
      }
      worst = larger(worst, {update_norm, cell.index});
      for (int s = 0; s < fields; s++) {
        previous_state[idx][s] = state[s];
      }
    }
  }
  if (first_report) {
    first_report = false;
    return;
  }
  for (int s = 0; s < fields; s++) {
    max_updates[s] = max_delta[s];
  }

  static bool header = false;
  if (!header) {
    printf("%7s %12s %12s %12s %12s %12s  %s\n", "step", "drho", "du", "dv", "de", "dnu_tilde",
           "maxcell");
    header = true;
  }
  printf("%7d", step);
  for (int s = 0; s < fields; s++) {
    printf(" %12.6e", max_updates[s]);
  }
  const auto& cell = cfd::cell_by_id(worst.cell);
  printf("  #%d(%.4f,%.4f)\n", worst.cell, cell.center.x, cell.center.y);
}

} // namespace convergence
