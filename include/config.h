#pragma once
#include <filesystem>

namespace config {
enum class FlowModel { laminar, sa };

struct Settings {
  FlowModel model = FlowModel::sa;
  std::filesystem::path mesh_path, log_path, field_path;
  int max_steps = 0, dump_interval = 0, convergence_interval = 0;
  double cfl = 0.0;
};
inline Settings settings;
inline bool uses_sa() {
  return settings.model == FlowModel::sa;
}
inline const char* model_name() {
  return uses_sa() ? "sa" : "laminar";
}
void load(const std::filesystem::path& path);
} // namespace config
