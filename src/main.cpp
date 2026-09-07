#include "config.h"
#include "parallel.h"
#include "solver.h"
#include <cstdio>
#include <exception>
#include <filesystem>
#include <omp.h>
#include <stdexcept>

int main(int argc, char** argv) {
  try {
    if (argc > 2)
      throw std::runtime_error("Usage: purinecfd [config.json]");
    const char* thread_policy = parallel::configure_threads();
    config::load(argc == 2 ? argv[1] : "config.json");
    const auto& options = config::settings;
    std::filesystem::create_directories(options.log_path.parent_path());
    std::filesystem::create_directories(options.field_path);
    if (!std::freopen(options.log_path.c_str(), "w", stdout))
      throw std::runtime_error("Cannot open log: " + options.log_path.string());
    std::printf("Steady %s | CFL=%.2f | max_steps=%d | OpenMP threads=%d | thread_policy=%s\n",
                config::uses_sa() ? "SA-RANS" : "laminar Navier-Stokes", options.cfl,
                options.max_steps, omp_get_max_threads(), thread_policy);
    solve_steady_flow();
    if (std::fflush(stdout) != 0)
      throw std::runtime_error("Failed to flush the run log");
    return 0;
  } catch (const std::exception& error) {
    std::fprintf(stderr, "Error: %s\n", error.what());
    return 1;
  }
}
