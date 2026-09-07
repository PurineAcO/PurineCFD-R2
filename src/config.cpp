#include "config.h"
#include "physics.h"
#include <cmath>
#include <fstream>
#include <limits>
#include <nlohmann/json.hpp>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Json = nlohmann::json;

void keys(const Json& object, std::initializer_list<const char*> allowed) {
  if (!object.is_object())
    throw std::runtime_error("Expected a configuration object");
  const std::set<std::string> names(allowed.begin(), allowed.end());
  for (const auto& item : object.items()) {
    if (!names.count(item.key()))
      throw std::runtime_error("Unknown configuration key: " + item.key());
  }
  for (const auto& name : names) {
    if (!object.contains(name))
      throw std::runtime_error("Missing configuration key: " + name);
  }
}

double positive(const Json& object, const char* key) {
  const auto& value = object.at(key);
  if (!value.is_number())
    throw std::runtime_error(std::string(key) + " must be a number");
  const double number = value.get<double>();
  if (!std::isfinite(number) || number <= 0.0)
    throw std::runtime_error(std::string(key) + " must be finite and positive");
  return number;
}

int count(const Json& object, const char* key) {
  const double value = positive(object, key);
  if (!object.at(key).is_number_integer() || value > std::numeric_limits<int>::max()) {
    throw std::runtime_error(std::string(key) + " must be a positive 32-bit integer");
  }
  return static_cast<int>(value);
}

std::filesystem::path path_value(const Json& object, const char* key,
                                 const std::filesystem::path& base) {
  const auto& value = object.at(key);
  if (!value.is_string())
    throw std::runtime_error(std::string(key) + " must be a path string");
  const auto text = value.get<std::string>();
  if (text.empty() || text.find('\0') != std::string::npos)
    throw std::runtime_error(std::string(key) + " is an invalid path");
  return std::filesystem::absolute(base / text).lexically_normal();
}
} // namespace

void config::load(const std::filesystem::path& path) {
  std::ifstream stream(path);
  if (!stream)
    throw std::runtime_error("Cannot open configuration: " + path.string());
  // Reject duplicate keys instead of allowing the last occurrence to silently win.
  std::vector<std::set<std::string>> object_keys;
  auto callback = [&](int depth, Json::parse_event_t event, Json& parsed) {
    if (depth > 16)
      throw std::runtime_error("Configuration nesting is too deep");
    if (event == Json::parse_event_t::object_start)
      object_keys.emplace_back();
    if (event == Json::parse_event_t::key &&
        !object_keys.back().insert(parsed.get<std::string>()).second) {
      throw std::runtime_error("Duplicate configuration key: " + parsed.get<std::string>());
    }
    if (event == Json::parse_event_t::object_end)
      object_keys.pop_back();
    return true;
  };
  const auto root = Json::parse(stream, callback);
  keys(root, {"io", "solver", "farfield"});
  const auto& io = root.at("io");
  const auto& solver = root.at("solver");
  const auto& far = root.at("farfield");
  keys(io, {"mesh", "log", "field"});
  keys(solver, {"max_steps", "cfl", "dump_interval", "convergence_interval"});
  keys(far, {"Ma", "T", "p"});

  Settings result;
  const auto base = std::filesystem::absolute(path).parent_path();
  result.mesh_path = path_value(io, "mesh", base);
  result.log_path = path_value(io, "log", base);
  result.field_path = path_value(io, "field", base);
  if (result.log_path == result.mesh_path ||
      result.log_path == std::filesystem::absolute(path).lexically_normal()) {
    throw std::runtime_error("The log path must not overwrite an input file");
  }
  result.max_steps = count(solver, "max_steps");
  result.dump_interval = count(solver, "dump_interval");
  result.convergence_interval = count(solver, "convergence_interval");
  if (result.convergence_interval < 2)
    throw std::runtime_error("convergence_interval must be at least 2");
  result.cfl = positive(solver, "cfl");
  const double ma = positive(far, "Ma");
  if (ma >= 1.0)
    throw std::runtime_error("This solver requires a subsonic farfield: 0 < Ma < 1");
  const double temperature = positive(far, "T");
  const double pressure = positive(far, "p");
  const double velocity = ma * sound_speed(temperature);
  const double density = pressure / (cfd::R * temperature);
  if (!std::isfinite(velocity) || !std::isfinite(density) || density <= 0.0) {
    throw std::runtime_error("Farfield values produce an invalid thermodynamic state");
  }
  settings = result;
  cfd::freestream = {velocity, 0.0, temperature, pressure};
}
