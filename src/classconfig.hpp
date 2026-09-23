#pragma once

#include <atomic>
#include <vector>
#include "config.hpp"
#include "physic.hpp"
#include <cmath>

namespace cc {

struct cell_class;      // 网格
struct face_class;      // 面

struct LSCBmatrix{
    double LU = 0;         // 上三角位置
    double RD = 0;         // 下三角位置
    double DC = 0;         // 对称位置
    double dxi[4] ={};     // delta x
    double dyi[4] ={};     // delta y
    double wi[4] = {};     // 范数权
};

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
    physics lowp;                  // 左插值物理量
    physics highp;                 // 右插值物理量

    // 结构化网格参数
    bool iswedir = false;           // 东西面指示
    cell_class *low,*high;          // 高低侧网格指针

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
    LSCBmatrix LSCB;             // LSCB梯度预处理矩阵

    // 用于结构化网格选项
    short east = -1,west = -1,north = -1,south = -1; // 东/西/南/北侧面在本格 faces 中的下标
    face_class* eastf = nullptr;  // 东侧邻接面
    face_class* westf = nullptr;  // 西侧邻接面
    face_class* northf = nullptr; // 北侧邻接面
    face_class* southf = nullptr; // 南侧邻接面
    int s = 0,n = 0;              // 环向/径向索引

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
inline std::vector<cell_class> GhostList;   // 虚网格
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

namespace cc {

inline node_class::node_class(int number_,double x_,double y_):number(number_),x(x_),y(y_){}

inline face_class::face_class(int index_,int p1_,int p2_,int c1_,int c2_,short type_)
    :index(index_),type(type_),cell_1(c1_),cell_2(c2_){
    node[0] = &NodeList[p1_-1];
    node[1] = &NodeList[p2_-1];
    mid = vec2{0.5*(node[0]->x + node[1]->x),0.5*(node[0]->y + node[1]->y)};
    nor = vec2{node[0]->y - node[1]->y,node[1]->x - node[0]->x};
    len = std::hypot(nor.x,nor.y);
}

inline cell_class::cell_class(int index_,int f1_,int f2_,int f3_,int f4_)
    :index(index_),face{f1_,f2_,f3_,f4_}{}

inline void face_class::face_physic_mid(){
    if(type != INTER){
        return;
    }
    phy.rho = 0.5*(nei[0]->phy.rho + nei[1]->phy.rho);
    phy.u   = 0.5*(nei[0]->phy.u   + nei[1]->phy.u);
    phy.v   = 0.5*(nei[0]->phy.v   + nei[1]->phy.v);
    phy.T   = 0.5*(nei[0]->phy.T   + nei[1]->phy.T);
    tur.miubl = 0.5*(nei[0]->tur.miubl + nei[1]->tur.miubl);
}

inline void face_class::form_physic(){
    phy.e = get_energy(phy);
    phy.p = cc::R*phy.rho*phy.T;
    phy.a = get_sonic_velocity(phy.T);
}

inline void cell_class::form_physic(){
    phy.e = get_energy(phy);
    phy.p = cc::R*phy.rho*phy.T;
    phy.a = get_sonic_velocity(phy.T);
}

inline void cell_class::form_conservative(){
    conser[0] = phy.rho;
    conser[1] = phy.rho*phy.u;
    conser[2] = phy.rho*phy.v;
    conser[3] = phy.rho*phy.e;
}

inline void cell_class::face_normal_out(){
    for(int i=0;i<ecnt;i++){
        fnorm[i] = dot(faces[i]->nor,
                       vec2{faces[i]->mid.x - center.x,faces[i]->mid.y - center.y}) > 0;
    }
}

inline void cell_class::reform(){
    phy.rho = conser[0];
    phy.u = conser[1]/conser[0];
    phy.v = conser[2]/conser[0];
    phy.e = conser[3]/conser[0];
    phy.T = (phy.e - 0.5*(phy.u*phy.u + phy.v*phy.v))/cc::Cv;
}

inline void cell_class::copyconver(){
    for(int i=0;i<4;i++){
        conserformer[i] = conser[i];
    }
    tur.miubl_former = tur.miubl;
}

inline cell_class& gotocell(int number){ return CellList[number-1]; }

inline face_class& gotoface(int number){ return FaceList[number-1]; }

inline face_class* link_face(int number){ return &FaceList[number-1]; }

inline cell_class* link_cell(int number){ return number <= 0 ? nullptr : &CellList[number-1]; }

inline cell_class* boundary_findcell(face_class* face){return face->nei[0] ? face->nei[0] : face->nei[1];}

inline bool field_ok(const cell_class& cell){
    const physics& phy = cell.phy;
    return std::isfinite(phy.rho) && phy.rho > 0.0 &&
           std::isfinite(phy.u) && std::isfinite(phy.v) &&
           std::isfinite(phy.T) && phy.T > 0.0 &&
           std::isfinite(phy.a) && std::isfinite(phy.p) && std::isfinite(phy.e) &&
           std::isfinite(cell.tur.miubl) && cell.tur.miubl >= 0.0;
}

inline void field_mark(const cell_class& cell){
    if(std::isfinite(cell.conser[0]) && cell.conser[0] > 0.0 &&
       std::isfinite(cell.conser[1]) && std::isfinite(cell.conser[2]) &&
       std::isfinite(cell.conser[3]) && std::isfinite(cell.tur.miubl)){
        return;
    }
    int current = field_bad_cell.load(std::memory_order_relaxed);
    while((current == 0 || cell.index < current) &&
          !field_bad_cell.compare_exchange_weak(current,cell.index,
                                                std::memory_order_relaxed)){
    }
}

}
