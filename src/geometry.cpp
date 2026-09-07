#include "geometry.h"
#include "mesh.h"
#include "mesh_reader.h"
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>

namespace {
bool contains_node(int num[], int number) {
  for (int i = 0; i < 4; i++) {
    if (number == num[i]) {
      return true;
    }
  }
  return false;
}

double point_distance(cfd::Vector2 a, cfd::Vector2 b) {
  return std::sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}

void collect_cell_nodes(cfd::Cell& cell) {
  short place = 0;
  for (int i = 0; i < cell.face_count; i++) {
    if (!contains_node(cell.node_indices, cell.faces[i]->nodes[0]->number - 1)) {
      cell.node_indices[place] = cell.faces[i]->nodes[0]->number - 1;
      place++;
    }
    if (!contains_node(cell.node_indices, cell.faces[i]->nodes[1]->number - 1)) {
      cell.node_indices[place] = cell.faces[i]->nodes[1]->number - 1;
      place++;
    }
    if (place == cell.face_count) {
      break;
    }
  }
}

double triangle_area(double x1, double y1, double x2, double y2, double x3, double y3) {
  return std::abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)) * 0.5;
}

cfd::Vector2 triangle_center(double x1, double y1, double x2, double y2, double x3, double y3) {
  return {(x1 + x2 + x3) / 3, (y1 + y2 + y3) / 3};
}

void compute_cell_area(cfd::Cell& cell) {
  const auto& p0 = cfd::nodes[cell.node_indices[0]];
  const auto& p1 = cfd::nodes[cell.node_indices[1]];
  const auto& p2 = cfd::nodes[cell.node_indices[2]];

  const auto& p3 = cfd::nodes[cell.node_indices[3]];
  // 四个顶点任选三个组成四个三角形，其面积和等于凸四边形面积的两倍。
  double S[4];
  S[0] = triangle_area(p0.x, p0.y, p1.x, p1.y, p2.x, p2.y);
  S[1] = triangle_area(p0.x, p0.y, p1.x, p1.y, p3.x, p3.y);
  S[2] = triangle_area(p0.x, p0.y, p2.x, p2.y, p3.x, p3.y);
  S[3] = triangle_area(p1.x, p1.y, p2.x, p2.y, p3.x, p3.y);
  cell.volume = (S[0] + S[1] + S[2] + S[3]) * 0.5;
}

void compute_cell_centroid(cfd::Cell& cell) {
  const auto& p0 = cfd::nodes[cell.node_indices[0]];
  const auto& p1 = cfd::nodes[cell.node_indices[1]];
  const auto& p2 = cfd::nodes[cell.node_indices[2]];

  const auto& p3 = cfd::nodes[cell.node_indices[3]];
  // 将四个三角形的形心按面积加权，得到凸四边形的形心。
  double S[4];
  cfd::Vector2 G[4];

  S[0] = triangle_area(p0.x, p0.y, p1.x, p1.y, p2.x, p2.y);
  G[0] = triangle_center(p0.x, p0.y, p1.x, p1.y, p2.x, p2.y);
  S[1] = triangle_area(p0.x, p0.y, p1.x, p1.y, p3.x, p3.y);
  G[1] = triangle_center(p0.x, p0.y, p1.x, p1.y, p3.x, p3.y);
  S[2] = triangle_area(p0.x, p0.y, p2.x, p2.y, p3.x, p3.y);
  G[2] = triangle_center(p0.x, p0.y, p2.x, p2.y, p3.x, p3.y);
  S[3] = triangle_area(p1.x, p1.y, p2.x, p2.y, p3.x, p3.y);
  G[3] = triangle_center(p1.x, p1.y, p2.x, p2.y, p3.x, p3.y);

  cell.center = {(S[0] * G[0].x + S[1] * G[1].x + S[2] * G[2].x + S[3] * G[3].x) /
                     (S[0] + S[1] + S[2] + S[3]),
                 (S[0] * G[0].y + S[1] * G[1].y + S[2] * G[2].y + S[3] * G[3].y) /
                     (S[0] + S[1] + S[2] + S[3])};
}

} // namespace

void cache_wall_distance(cfd::Cell& cell) {
  double distance = std::numeric_limits<double>::infinity();
  for (cfd::Face* wall : cfd::wall_faces) {
    double length = point_distance(cell.center, wall->midpoint);
    distance = std::min(distance, length);
  }
  cell.turbulence.inverse_wall_distance_squared = 1.0 / (distance * distance);
}

void initialize_geometry() {
  connect_mesh();
  for (cfd::Cell& cell : cfd::cells) {
    collect_cell_nodes(cell);
    compute_cell_area(cell);
    compute_cell_centroid(cell);
    cell.orient_face_normals();
    if (!std::isfinite(cell.volume) || cell.volume <= 0.0) {
      throw std::runtime_error("Invalid cell area");
    }
    cell.inverse_volume = 1.0 / cell.volume;
    for (const auto* face : cell.faces) {
      cell.projected_face_sum.x += 0.5 * std::abs(face->area_normal.x);
      cell.projected_face_sum.y += 0.5 * std::abs(face->area_normal.y);
    }
  }
}
