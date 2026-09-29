#pragma once

#include <atomic>
#include <vector>
#include "config.hpp"
#include <cmath>

#define allface(cell) for(int i=0;i<cell.ecnt;i++)
#define allcell for(int i=0;i<cc::cell_num;i++)
#define icell(i) cc::CellList[i]
#define allfac for(int i=0;i<cc::face_num;i++)
#define iface(i) cc::FaceList[i]

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
    int number = 0;                     // 节点编号
    double x = 0.0, y = 0.0;            // 节点坐标
    node_class() = default;
    node_class(int number_,double x_,double y_);
};

struct face_class{
    int index = 0;                  // 面编号
    short type = INTER;             // 面类型
    node_class* node[2] = {};      // 面邻接点指针
    vec2 mid = {0.0,0.0};   // 面中点坐标
    vec2 nor = {0.0,0.0};   // 带面长法向量
    bool outer;                    // 面法向方向01指示
    double len = 0.0;               // 面长度
    int cell_1 = -1, cell_2 = -1;  // 面邻接网格编号
    cell_class* nei[2] = {};       // 面邻接网格指针
    vecp phy;                       // 物理量
    vecgrad phgrad;                 // 面法向梯度
    struct otphy otphy;             // 引申物理量
    double un;                     // 面法向速度
    turbulence tur;                // 湍流
    vecp phynei[2];                // 左右插值物理量
    struct otphy otnei[2];         // 左右引申物理量

    // 结构化网格参数
    // bool iswedir = false;       // 面方向指示
    vec4 convect;               // 无粘通量
    vec4 visflux;               // 黏性通量

    face_class() = default;
    face_class(int index_,int p1_,int p2_,int c1_,int c2_,short type_);

    // 法向量单位化
    vec2 length1_nor();
    // 指示面法向方向
    void normal_out();
    // 面上中心差分插值
    void face_physic_mid();
    // 根据基本物理量形成引申物理量
    void form_otherphy();
    // 通用级jacobi变换,已经带了面长度
    double toface_jacobi(double F,double G) const;
    cc::vec4 toface_jacobi(const cc::vec4 &F,const cc::vec4 &G) const;
    template<int N> void toface_jacobi(const double (&F)[N],const double (&G)[N],double (&out)[N]) const;
    template<int N1,int N2> void toface_jacobi(const double (&F)[N1][N2],const double (&G)[N1][N2],double (&out)[N1][N2]) const;
};

struct cell_class{
    int index = 0;                      // 编号
    int ecnt = 4;                       // 面邻接边个数
    int face[4] = {};                   // 邻接面编号
    int node[4] = {-1,-1,-1,-1}; // 邻接点下标
    bool fnorm[4] = {};                 // 邻接面外法向标记
    cell_class* nei[4] = {};            // 邻接网格指针
    face_class* faces[4] = {};          // 邻接面指针
    double vol = 0.0;                   // 体积
    vec2 center;                        // 中心坐标

    vecp phy;                    // 物理量
    struct otphy otphy;                 // 引申物理量
    vec4 conser;                 // 守恒量
    vec4 conserformer;           // 前期守恒量
    vec4 lastconser;             // 上一时间步守恒量
    vec4 convect;                // 无粘对流项
    vec4 visflux;                // 黏性通量
    vecgrad phgrad;              // 物理量梯度
    dissipation diss;            // 耗散项
    double localdt = 0.0;        // 当地时间步长
    turbulence tur;              // 湍流
    LSCBmatrix LSCB;             // LSCB梯度预处理矩阵

    // 东/西/南/北侧面在本格 faces 中的下标
    short east = -1,west = -1,north = -1,south = -1; 
    face_class* eastf = nullptr;  // 东侧邻接面
    face_class* westf = nullptr;  // 西侧邻接面
    face_class* northf = nullptr; // 北侧邻接面
    face_class* southf = nullptr; // 南侧邻接面
    int s = 0,n = 0;              // 环向/径向索引

    cell_class() = default;
    cell_class(int index_,int f1_,int f2_,int f3_,int f4_);

    // 由ρ,u,v,T形成e,p,a
    void form_otherphy();
    // 形成守恒量
    void form_conservative();
    // 找到邻接面外法向
    void face_normal_out();
    // 由守恒量恢复ρ,u,v,T
    void prim();
    // 保存本步RK的基准状态
    void copyconver();
    // 保存本时间步的基准状态
    void copyconver_time();
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

inline void face_class::form_otherphy(){
    otphy.form_otphy(phy);
    un = toface_jacobi(phy.u,phy.v)/nor.norm();
}

inline void cell_class::form_otherphy(){
    otphy.form_otphy(phy);
}

inline void cell_class::form_conservative(){
    conser = vec4(phy.rho,phy.rho*phy.u,phy.rho*phy.v,phy.rho*otphy.e);
}

inline void cell_class::face_normal_out(){
    for(int i=0;i<ecnt;i++){
        fnorm[i] = dot(faces[i]->nor,faces[i]->mid-center) > 0;
    }
}

inline void face_class::normal_out(){
    outer = dot(nor,nei[1]->center-nei[0]->center) > 0;
}

inline void cell_class::prim(){
    phy.rho = conser.c;
    phy.u = conser.x/conser.c;
    phy.v = conser.y/conser.c;
    double e = conser.e/conser.c;
    phy.T = (e - 0.5*(phy.u*phy.u + phy.v*phy.v))/cc::Cv;
    form_otherphy();
}

inline void cell_class::copyconver(){
    conserformer = conser;
    if(!cc::urans)tur.miubl_former = tur.miubl;
}

inline void cell_class::copyconver_time(){
    lastconser = conser;
    if(cc::urans)tur.miubl_former = tur.miubl;
}

inline cell_class& gotocell(int number){ return CellList[number-1]; }

inline face_class& gotoface(int number){ return FaceList[number-1]; }

inline face_class* link_face(int number){ return &FaceList[number-1]; }

inline cell_class* link_cell(int number){ return number <= 0 ? nullptr : &CellList[number-1]; }

inline cell_class* boundary_findcell(face_class* face){return face->nei[0] ? face->nei[0] : face->nei[1];}

inline bool field_ok(const cell_class& cell){
    const vecp& phy = cell.phy;
    return std::isfinite(phy.rho) && phy.rho > 0.0 &&
           std::isfinite(phy.u) && std::isfinite(phy.v) &&
           std::isfinite(phy.T) && phy.T > 0.0 &&
           std::isfinite(cell.otphy.a) && std::isfinite(cell.otphy.p) && std::isfinite(cell.otphy.e) &&
           std::isfinite(cell.tur.miubl) && cell.tur.miubl >= 0.0;
}

inline void field_mark(const cell_class& cell){
    if(std::isfinite(cell.conser.c) && cell.conser.c > 0.0 &&
       std::isfinite(cell.conser.x) && std::isfinite(cell.conser.y) &&
       std::isfinite(cell.conser.e) && std::isfinite(cell.tur.miubl)){
        return;
    }
    int current = field_bad_cell.load(std::memory_order_relaxed);
    while((current == 0 || cell.index < current) &&
          !field_bad_cell.compare_exchange_weak(current,cell.index,
                                                std::memory_order_relaxed)){
    }
}

inline double face_class::toface_jacobi(double F,double G)const{
    return nor.x * F + nor.y * G;
}

inline cc::vec4 face_class::toface_jacobi(const cc::vec4 &F,const cc::vec4 &G) const{
    return nor.x * F + nor.y * G;
}

template<int N> void face_class::toface_jacobi(const double (&F)[N],const double (&G)[N],double (&out)[N]) const{
    for(int i=0;i<N;i++) out[i] = nor.x * F[i] + nor.y * G[i];
}

template<int N1,int N2> 
void face_class::toface_jacobi(const double (&F)[N1][N2],const double (&G)[N1][N2],double (&out)[N1][N2]) const{
    for(int i=0;i<N1;i++) for(int j=0;j<N2;j++) out[i][j] = nor.x*F[i][j] + nor.y * G[i][j];
}

inline vec2 face_class::length1_nor(){
    return vec2{nor.x/len,nor.y/len};
}


}
