#include "parallel.h"
#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <omp.h>
#include <sched.h>
#include <set>
#include <stdexcept>
#include <string>
#include <system_error>
#include <unistd.h>
#include <utility>

namespace {
int physical_cores() {
  const long configured = sysconf(_SC_NPROCESSORS_CONF);
  if (configured <= 0)
    throw std::runtime_error("Cannot determine the configured CPU count");
  const int capacity = static_cast<int>(configured);
  const auto release = [](cpu_set_t* mask) { CPU_FREE(mask); };
  const std::unique_ptr<cpu_set_t, decltype(release)> mask(CPU_ALLOC(capacity), release);
  const auto bytes = CPU_ALLOC_SIZE(capacity);
  if (!mask)
    throw std::bad_alloc();
  if (sched_getaffinity(0, bytes, mask.get()) != 0)
    throw std::system_error(errno, std::generic_category(), "Cannot read CPU affinity");

  std::set<std::pair<int, int>> cores;
  for (int cpu = 0; cpu < capacity; ++cpu) {
    if (!CPU_ISSET_S(cpu, bytes, mask.get()))
      continue;
    const auto base = "/sys/devices/system/cpu/cpu" + std::to_string(cpu) + "/topology/";
    int socket = -1, core = -1;
    std::ifstream socket_file(base + "physical_package_id");
    std::ifstream core_file(base + "core_id");
    if (!(socket_file >> socket) || !(core_file >> core) || socket < 0 || core < 0) {
      throw std::runtime_error("Cannot read physical topology for CPU " + std::to_string(cpu));
    }
    cores.emplace(socket, core);
  }
  if (cores.empty())
    throw std::runtime_error("No physical CPU cores available");
  return static_cast<int>(cores.size());
}
} // namespace

const char* parallel::configure_threads() {
  omp_set_dynamic(0);
  if (const char* requested = std::getenv("OMP_NUM_THREADS"); requested && *requested)
    return "OMP_NUM_THREADS";
  omp_set_num_threads(std::min(physical_cores(), omp_get_thread_limit()));
  return "available-physical-cores";
}
