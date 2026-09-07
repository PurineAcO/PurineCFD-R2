#pragma once
#include "flow.h"
#include <vector>

namespace cfd {

// 文件中的节点、面、单元编号从 1 开始；vector 下标从 0 开始。

inline int cell_count = 0;
inline int face_count = 0;
enum class BoundaryType { interior, wall, farfield };

struct Cell;
struct Face;

struct Node {
  int number = 0;          // 节点编号
  double x = 0.0, y = 0.0; // 节点坐标
  Node(int number_, double x_, double y_);
};

struct Face {
  int index = 0;                               // 面编号
  BoundaryType type = BoundaryType::interior;  // 面类型
  Node* nodes[2] = {};                         // 面两端的节点
  Vector2 midpoint = {0.0, 0.0};               // 面中点坐标
  double length = 0.0;                         // 面长，m；初始化后固定
  Vector2 area_normal = {0.0, 0.0};            // n*Δs，含面长，方向由节点顺序确定
  int first_cell_id = -1, second_cell_id = -1; // 面两侧的单元编号
  Cell* cells[2] = {};                         // 面两侧的单元；边界外侧为 nullptr
  FlowState flow;                              // 流动状态
  TurbulenceState turbulence;                  // SA 模型变量

  // 每个 RK 阶段在面变量及梯度更新后计算；每个面只写一次。
  double volume_flux = 0.0; // 单位厚度的体积流量 (u·n)*Δs，m²/s
  double spectral_radius = 0.0;
  double sa_diffusivity = 0.0;
  double sa_gradient_flux = 0.0; // (∇ν̃·n)*Δs
  double convective_flux[4] = {};
  double viscous_flux[4] = {};

  Face(int index_, int p1_, int p2_, int c1_, int c2_, BoundaryType type_);

  // 内部面原始变量取两侧单元的算术平均；边界面由边界条件赋值
  void interpolate_primitives();
  // 由面上的 ρ、u、v、T 计算总能量 E、压力 p 和声速 a。
  void update_thermodynamics();
};

struct Cell {
  int index = 0;                          // 编号
  static constexpr int face_count = 4;    // 四边形单元的面数
  int face_ids[4] = {};                   // 邻接面编号
  int node_indices[4] = {-1, -1, -1, -1}; // 节点 vector 下标，从 0 开始
  bool normal_points_outward[4] = {};     // true 表示该面的 area_normal 指向本单元外侧
  Cell* neighbors[4] = {};                // 相邻单元；边界处为 nullptr
  Face* faces[4] = {};                    // 邻接面指针
  double inverse_volume = 0.0;
  Vector2 projected_face_sum;
  double volume = 0.0; // 二维单元面积（按单位厚度计），m²
  Vector2 center;      // 中心坐标

  FlowState flow;                       // 流动状态
  double conservative[4] = {};          // Q = [ρ, ρu, ρv, ρE]
  double previous_conservative[4] = {}; // Qⁿ：本步 RK 开始时保存的状态
  double viscous_flux[4] = {};
  double convective_flux[4] = {};   // 无粘对流项
  JstDissipation dissipation_terms; // 耗散项
  double local_dt = 0.0;            // 局部伪时间步长，s
  TurbulenceState turbulence;       // SA 模型变量

  Cell(int index_, int f1_, int f2_, int f3_, int f4_);

  // 根据原始变量计算 E、p、a
  void update_thermodynamics();
  // 由 ρ、u、v、E 计算 Q = [ρ, ρu, ρv, ρE]。
  void update_conservative();
  // 根据单元中心与面中点的位置，判断面法向是否朝外。
  void orient_face_normals();
  // 由 Q 恢复 ρ、u、v、T、E
  void recover_primitives();
  // 保存本步 RK 的基准状态 Qⁿ
  void save_previous_state();
};

inline std::vector<Node> nodes;           // 节点
inline std::vector<Cell> cells;           // 网格
inline std::vector<Face> faces;           // 面
inline std::vector<Face*> wall_faces;     // 壁面
inline std::vector<Face*> farfield_faces; // 远场

// 按文件编号访问单元，编号越界时抛出异常。
Cell& cell_by_id(int number);
// 按文件编号访问面，编号越界时抛出异常。
Face& face_by_id(int number);
// 仅用于边界面：返回唯一的内部单元
Cell* boundary_cell(Face* face);

} // namespace cfd
