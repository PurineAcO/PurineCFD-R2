#include "mesh.h"
#include "physics.h"
#include <cmath>
#include <stdexcept>
#include <string>

namespace cfd {

Node::Node(int number_, double x_, double y_) : number(number_), x(x_), y(y_) {}

Face::Face(int index_, int p1_, int p2_, int c1_, int c2_, BoundaryType type_)
    : index(index_), type(type_), first_cell_id(c1_), second_cell_id(c2_) {
  nodes[0] = &cfd::nodes[p1_ - 1];
  nodes[1] = &cfd::nodes[p2_ - 1];
  midpoint = {0.5 * (cfd::nodes[p1_ - 1].x + cfd::nodes[p2_ - 1].x),
              0.5 * (cfd::nodes[p1_ - 1].y + cfd::nodes[p2_ - 1].y)};
  area_normal = {cfd::nodes[p1_ - 1].y - cfd::nodes[p2_ - 1].y,
                 cfd::nodes[p2_ - 1].x - cfd::nodes[p1_ - 1].x};
  length = std::hypot(area_normal.x, area_normal.y);
}

Cell::Cell(int index_, int f1_, int f2_, int f3_, int f4_)
    : index(index_), face_ids{f1_, f2_, f3_, f4_} {}

void Cell::update_conservative() {
  conservative[0] = flow.rho;
  conservative[1] = flow.rho * flow.u;
  conservative[2] = flow.rho * flow.v;
  conservative[3] = flow.rho * flow.e;
}

void Cell::orient_face_normals() {
  for (int i = 0; i < face_count; i++) {
    normal_points_outward[i] = dot(faces[i]->area_normal, {faces[i]->midpoint.x - center.x,
                                                           faces[i]->midpoint.y - center.y}) > 0;
  }
}

void Cell::update_thermodynamics() {
  flow.e = cfd::Cv * flow.T + 0.5 * (flow.u * flow.u + flow.v * flow.v);
  flow.p = cfd::R * flow.rho * flow.T;
  flow.a = sound_speed(flow.T);
}

void Cell::recover_primitives() {
  flow.rho = conservative[0];
  flow.u = conservative[1] / conservative[0];
  flow.v = conservative[2] / conservative[0];
  flow.e = conservative[3] / conservative[0];
  flow.T = (flow.e - 0.5 * (flow.u * flow.u + flow.v * flow.v)) / cfd::Cv;
}

void Cell::save_previous_state() {
  for (int i = 0; i < 4; i++) {
    previous_conservative[i] = conservative[i];
  }
}

void Face::interpolate_primitives() {
  if (type == BoundaryType::interior) {
    flow.rho = 0.5 * (cells[0]->flow.rho + cells[1]->flow.rho);
    flow.u = 0.5 * (cells[0]->flow.u + cells[1]->flow.u);
    flow.v = 0.5 * (cells[0]->flow.v + cells[1]->flow.v);
    flow.T = 0.5 * (cells[0]->flow.T + cells[1]->flow.T);
  }
}

void Face::update_thermodynamics() {
  flow.e = cfd::Cv * flow.T + 0.5 * (flow.u * flow.u + flow.v * flow.v);
  flow.p = cfd::R * flow.rho * flow.T;
  flow.a = sound_speed(flow.T);
}

Cell& cell_by_id(int number) {
  if (number < 1 || number > cfd::cell_count) {
    throw std::out_of_range("Cell index out of range: " + std::to_string(number));
  }
  return cfd::cells[number - 1];
}

Face& face_by_id(int number) {
  if (number < 1 || number > cfd::face_count) {
    throw std::out_of_range("Face index out of range: " + std::to_string(number));
  }
  return cfd::faces[number - 1];
}

Cell* boundary_cell(Face* face) {
  return face->cells[0] ? face->cells[0] : face->cells[1];
}

} // namespace cfd
