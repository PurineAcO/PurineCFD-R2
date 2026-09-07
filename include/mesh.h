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
  int number = 0;                          // 节点编号
  double x = 0.0, y = 0.0;                 // 节点坐标
  Node(int number_, double x_, double y_); // 节点构造器
};

struct Face {
  int index = 0;                               // 面编号
  BoundaryType type = BoundaryType::interior;  // 面类型
  Node* nodes[2] = {};                         // 面邻接点指针
  Vector2 midpoint = {0.0, 0.0};               // 面中点坐标
  double length = 0.0;                         // 面长，m；初始化后固定
  Vector2 area_normal = {0.0, 0.0};            // n*Δs，含面长，方向由节点顺序确定
  int first_cell_id = -1, second_cell_id = -1; // 面邻接网格编号
  Cell* cells[2] = {};                         // 面邻接网格指针
  FlowState flow;                              // 物理量
  TurbulenceState turbulence;                  // 湍流

  // 每个 RK 阶段在面变量及梯度更新后计算；每个面只写一次。
  double volume_flux = 0.0; // (u·n)*Δs，不是单位法向速度
  double spectral_radius = 0.0;
  double sa_diffusivity = 0.0;
  double sa_gradient_flux = 0.0; // (∇ν̃·n)*Δs
  double convective_flux[4] = {};
  double viscous_flux[4] = {};

  Face(int index_, int p1_, int p2_, int c1_, int c2_, BoundaryType type_); // 面构造

  // 内部面原始变量取两侧单元的算术平均；边界面由边界条件赋值
  void interpolate_primitives();
  // 形成面上所有物理量
  void update_thermodynamics();
};

struct Cell {
  int index = 0;                          // 编号
  static constexpr int face_count = 4;    // 面邻接边个数
  int face_ids[4] = {};                   // 邻接面编号
  int node_indices[4] = {-1, -1, -1, -1}; // 节点 vector 下标，从 0 开始
  bool normal_points_outward[4] = {};     // 邻接面外法向标记
  Cell* neighbors[4] = {};                // 相邻单元；边界处为 nullptr
  Face* faces[4] = {};                    // 邻接面指针
  double inverse_volume = 0.0;
  Vector2 projected_face_sum;
  double volume = 0.0; // 二维单元面积（按单位厚度计），m²
  Vector2 center;      // 中心坐标

  FlowState flow;                       // 物理量
  double conservative[4] = {};          // Q = [ρ, ρu, ρv, ρE]
  double previous_conservative[4] = {}; // Qⁿ：本步 RK 开始时保存的状态
  double viscous_flux[4] = {};
  double convective_flux[4] = {};   // 无粘对流项
  JstDissipation dissipation_terms; // 耗散项
  double local_dt = 0.0;            // 当地时间步长
  TurbulenceState turbulence;       // 湍流

  Cell(int index_, int f1_, int f2_, int f3_, int f4_); // 网格构造器

  // 根据原始变量计算 E、p、a
  void update_thermodynamics();
  // 形成守恒量
  void update_conservative();
  // 找到邻接面法向
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

// 安全访问网格
Cell& cell_by_id(int number);
// 安全访问邻接边
Face& face_by_id(int number);
// 仅用于边界面：返回唯一的内部单元
Cell* boundary_cell(Face* face);

} // namespace cfd
