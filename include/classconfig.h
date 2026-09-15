#pragma once
#include <atomic>
#include <vector>
#include "config.h"

namespace cc {

struct cell_class;      // 网格
struct face_class;      // 面

struct node_class{
    int number = 0;         // 节点编号
    double x = 0.0, y = 0.0;// 节点坐标
    node_class() = default;
    node_class(int number_,double x_,double y_);// 节点构造器
};

struct face_class{
    int index = 0;      // 面编号
    short type = INTER; // 面类型
    node_class* node[2] = {};      // 面邻接点指针
    vec2 mid = {0.0,0.0};   // 面中点坐标
    vec2 nor = {0.0,0.0};   // 面法向*面长
    double len = 0.0;       // 面长度
    int cell_1 = -1, cell_2 = -1;  // 面邻接网格编号
    cell_class* nei[2] = {};       // 面邻接网格指针
    physics phy;                   // 物理量
    turbulence tur;                // 湍流

    double volflux = 0.0;    // 单位厚度体积流量 (u·n)*Δs
    double lam = 0.0;        // 谱半径 |u·n|*Δs + a*Δs
    double coef = 0.0;       // ν̃的扩散系数
    double turflux = 0.0;    // (∇ν̃·n)*Δs
    double convect[4] = {};  // 无粘通量
    double visflux[4] = {};  // 黏性通量

    face_class() = default;
    face_class(int index_,int p1_,int p2_,int c1_,int c2_,short type_);// 面构造

    // 面上中心差分插值
    void face_physic_mid();
    // 面值由ρ,u,v,T形成e,p,a
    void form_physic();
};

struct cell_class{
    int index = 0;     // 编号
    int ecnt = 4;      // 面邻接边个数
    int face[4] = {};  // 邻接面编号
    int node[4] = {-1,-1,-1,-1}; // 邻接点下标
    bool fnorm[4] = {}; // 邻接面外法向标记
    cell_class* nei[4] = {};   // 邻接网格指针
    face_class* faces[4] = {}; // 邻接面指针
    double vol = 0.0;   // 体积(二维按单位厚度计)
    double invvol = 0.0;// 1/vol
    vec2 center;        // 中心坐标
    vec2 proj;          // 各方向投影面积和,用于当地时间步长

    physics phy;        // 物理量
    double conser[4] = {};       // 守恒量
    double conserformer[4] = {}; // 前期守恒量
    double convect[4] = {};      // 无粘对流项
    double visflux[4] = {};      // 黏性通量
    dissipation diss;            // 耗散项
    double localdt = 0.0;        // 当地时间步长
    turbulence tur;              // 湍流

    cell_class() = default;
    cell_class(int index_,int f1_,int f2_,int f3_,int f4_);// 网格构造器

    // 由ρ,u,v,T形成e,p,a
    void form_physic();
    // 形成守恒量
    void form_conservative();
    // 找到邻接面外法向
    void face_normal_out();
    // 由守恒量恢复ρ,u,v,T
    void reform();
    // 保存本步RK的基准状态
    void copyconver();
};

inline std::vector<node_class> NodeList;    // 节点
inline std::vector<cell_class> CellList;    // 网格
inline std::vector<face_class> FaceList;    // 面
inline std::vector<face_class*> WallFaces;  // 壁面
inline std::vector<face_class*> FarFaces;   // 远场

inline std::atomic<int> field_bad_cell{0};  // 首个异常网格编号,0表示正常

// 按文件编号访问网格
cell_class& gotocell(int number);
// 按文件编号访问面
face_class& gotoface(int number);
// 链接到面
face_class* link_face(int number);
// 链接到网格
cell_class* link_cell(int number);
// 对于边界面找到内部网格
cell_class* boundary_findcell(face_class* face);
// 物理量是否有效
bool field_ok(const cell_class& cell);
// 记录异常网格编号(取最小), 用于并行循环中标记
void field_mark(const cell_class& cell);

}
