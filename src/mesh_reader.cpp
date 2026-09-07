#include "mesh_reader.h"
#include "mesh.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
class MeshReader {
public:
  explicit MeshReader(const std::filesystem::path& path) : input(path) {
    if (!input)
      throw std::runtime_error("Cannot open mesh: " + path.string());
  }
  std::string line() {
    std::string text;
    if (!std::getline(input, text))
      fail("Unexpected end of mesh");
    ++line_number;
    if (!text.empty() && text.back() == '\r')
      text.pop_back();
    return text;
  }
  void expect(const std::string& expected) {
    if (line() != expected)
      fail("Expected " + expected);
  }
  template <typename... Values> void parse(const std::string& text, Values&... values) {
    std::istringstream row(text);
    if (!(row >> ... >> values))
      fail("Invalid mesh row");
    row >> std::ws;
    if (!row.eof())
      fail("Unexpected data after mesh row");
  }
  [[noreturn]] void fail(const std::string& message) const {
    throw std::runtime_error("Mesh line " + std::to_string(line_number) + ": " + message);
  }
  void finish() {
    input >> std::ws;
    if (!input.eof())
      fail("Unexpected data after cell section");
  }

private:
  std::ifstream input;
  int line_number = 0;
};
} // namespace

void read_mesh(const std::filesystem::path& path) {
  MeshReader reader(path);
  int node_count = 0, type_count = 0;
  reader.parse(reader.line(), node_count, cfd::face_count, cfd::cell_count, type_count);
  if (node_count < 4 || cfd::face_count < 4 || cfd::cell_count < 1 || type_count != 3)
    reader.fail("Expected positive counts and exactly INTER/WALL/FAR groups");
  std::map<std::string, cfd::BoundaryType> supported = {{"INTER", cfd::BoundaryType::interior},
                                                        {"WALL", cfd::BoundaryType::wall},
                                                        {"FAR", cfd::BoundaryType::farfield}};
  std::vector<std::pair<std::string, cfd::BoundaryType>> groups;
  std::set<std::string> names, kinds;
  for (int i = 0; i < type_count; ++i) {
    std::string name, equal, kind;
    reader.parse(reader.line(), name, equal, kind);
    if (equal != "=" || !supported.count(kind) || !names.insert(name).second ||
        !kinds.insert(kind).second)
      reader.fail("Invalid or duplicate boundary group");
    groups.emplace_back(name, supported.at(kind));
  }
  reader.expect("(node)");
  cfd::nodes.reserve(node_count);
  for (int i = 1; i <= node_count; ++i) {
    int index = 0;
    double x = 0.0, y = 0.0;
    reader.parse(reader.line(), index, x, y);
    if (index != i || !std::isfinite(x) || !std::isfinite(y))
      reader.fail("Invalid node index or coordinates");
    cfd::nodes.emplace_back(index, x, y);
  }
  reader.expect("(end)");
  reader.expect("(edge)");
  cfd::faces.reserve(cfd::face_count);
  std::set<int> face_ids;
  for (const auto& [name, type] : groups) {
    reader.expect(name);
    for (auto text = reader.line(); text != "(end)"; text = reader.line()) {
      int index = 0, a = 0, b = 0, left = 0, right = 0;
      reader.parse(text, index, a, b, left, right);
      if (index < 1 || index > cfd::face_count || !face_ids.insert(index).second)
        reader.fail("Invalid or duplicate face index");
      if (a < 1 || a > node_count || b < 1 || b > node_count || a == b)
        reader.fail("Invalid face nodes");
      if (left < 0 || left > cfd::cell_count || right < 0 || right > cfd::cell_count ||
          left == right)
        reader.fail("Invalid adjacent cell indices");
      const bool interior = left != 0 && right != 0;
      if (interior != (type == cfd::BoundaryType::interior))
        reader.fail("Boundary type disagrees with adjacency");
      const auto& pa = cfd::nodes[a - 1];
      const auto& pb = cfd::nodes[b - 1];
      if (!(std::hypot(pa.x - pb.x, pa.y - pb.y) > 0.0))
        reader.fail("Zero-length face");
      cfd::faces.emplace_back(index, a, b, left, right, type);
    }
  }
  if (cfd::faces.size() != static_cast<std::size_t>(cfd::face_count))
    reader.fail("Face count mismatch");
  std::sort(cfd::faces.begin(), cfd::faces.end(),
            [](const auto& a, const auto& b) { return a.index < b.index; });
  reader.expect("(end)");
  reader.expect("(cell)");
  cfd::cells.reserve(cfd::cell_count);
  std::vector<int> uses(cfd::face_count, 0);
  for (int i = 1; i <= cfd::cell_count; ++i) {
    int index = 0;
    std::array<int, 4> faces{};
    reader.parse(reader.line(), index, faces[0], faces[1], faces[2], faces[3]);
    if (index != i || std::set<int>(faces.begin(), faces.end()).size() != 4)
      reader.fail("Expected sequential quadrilateral cells with four distinct faces");
    std::map<int, int> degree;
    for (int id : faces) {
      if (id < 1 || id > cfd::face_count)
        reader.fail("Cell face index out of range");
      const auto& face = cfd::faces[id - 1];
      if (face.first_cell_id != i && face.second_cell_id != i)
        reader.fail("Cell/face adjacency mismatch");
      ++degree[face.nodes[0]->number];
      ++degree[face.nodes[1]->number];
      ++uses[id - 1];
    }
    if (degree.size() != 4 || std::any_of(degree.begin(), degree.end(),
                                          [](const auto& entry) { return entry.second != 2; }))
      reader.fail("Cell edges do not form a quadrilateral");
    // 逐边检查其余顶点是否位于同一侧，保证后续面积和形心公式适用于凸单元。
    for (int id : faces) {
      const auto& edge = cfd::faces[id - 1];
      const auto& a = *edge.nodes[0];
      const auto& b = *edge.nodes[1];
      double side = 0.0;
      for (const auto& [node_id, unused] : degree) {
        (void)unused;
        if (node_id == a.number || node_id == b.number)
          continue;
        const auto& c = cfd::nodes[node_id - 1];
        const double cross = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
        if (!std::isfinite(cross) || cross == 0.0 ||
            (side != 0.0 && std::signbit(side) != std::signbit(cross)))
          reader.fail("Cell is degenerate, concave or self-intersecting");
        side = cross;
      }
    }
    cfd::cells.emplace_back(index, faces[0], faces[1], faces[2], faces[3]);
  }
  for (const auto& face : cfd::faces) {
    const int expected = face.type == cfd::BoundaryType::interior ? 2 : 1;
    if (uses[face.index - 1] != expected)
      reader.fail("Face reference count mismatch");
  }
  reader.finish();
  std::printf("Mesh: nodes=%d faces=%d cells=%d\n", node_count, cfd::face_count, cfd::cell_count);
}

void connect_mesh() {
  for (auto& cell : cfd::cells) {
    for (int i = 0; i < cell.face_count; ++i)
      cell.faces[i] = &cfd::face_by_id(cell.face_ids[i]);
  }
  for (auto& face : cfd::faces) {
    face.cells[0] = face.first_cell_id == 0 ? nullptr : &cfd::cell_by_id(face.first_cell_id);
    face.cells[1] = face.second_cell_id == 0 ? nullptr : &cfd::cell_by_id(face.second_cell_id);
    if (face.type == cfd::BoundaryType::wall)
      cfd::wall_faces.push_back(&face);
    if (face.type == cfd::BoundaryType::farfield)
      cfd::farfield_faces.push_back(&face);
  }
  for (auto& cell : cfd::cells) {
    for (int i = 0; i < cell.face_count; ++i) {
      const auto& face = *cell.faces[i];
      cell.neighbors[i] = face.cells[0] == &cell ? face.cells[1] : face.cells[0];
    }
  }
  if (cfd::wall_faces.empty() || cfd::farfield_faces.empty())
    throw std::runtime_error("Mesh requires nonempty wall and farfield boundaries");
  std::printf("Boundaries: wall=%zu farfield=%zu\n", cfd::wall_faces.size(),
              cfd::farfield_faces.size());
}
