#pragma once
#include <filesystem>

namespace config {
struct Settings {
  std::filesystem::path mesh_path, log_path, field_path;
  int max_steps = 0, dump_interval = 0, convergence_interval = 0;
  double cfl = 0.0;
};
inline Settings settings;
void load(const std::filesystem::path& path);
} // namespace config
