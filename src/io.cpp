#include "io.h"
#include "config.h"
#include "mesh.h"
#include <cmath>
#include <cstdio>
#include <memory>
#include <stdexcept>

void dump_field(int step) {
  char filename[32];
  std::snprintf(filename, sizeof(filename), "step_%06d.dat", step);
  const auto path = config::settings.field_path / filename;
  using File = std::unique_ptr<std::FILE, decltype(&std::fclose)>;
  File file(std::fopen(path.c_str(), "w"), &std::fclose);
  if (!file)
    throw std::runtime_error("Cannot create field file: " + path.string());
  std::fprintf(file.get(), "TITLE=\"step %d model=%s\"\n", step, config::model_name());
  std::fprintf(file.get(),
               "VARIABLES=\"x\",\"y\",\"rho\",\"u\",\"v\",\"T\",\"p\",\"Ma\",\"nu_tilde\"\n");
  for (const auto& cell : cfd::cells) {
    std::fprintf(file.get(), "%.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e\n", cell.center.x,
                 cell.center.y, cell.flow.rho, cell.flow.u, cell.flow.v, cell.flow.T, cell.flow.p,
                 std::hypot(cell.flow.u, cell.flow.v) / cell.flow.a, cell.turbulence.nu_tilde);
  }
  if (std::ferror(file.get()))
    throw std::runtime_error("Failed to write field: " + path.string());
  if (std::fclose(file.release()) != 0)
    throw std::runtime_error("Failed to close field: " + path.string());
  std::printf("Field output: step %d %s\n", step, path.c_str());
}
