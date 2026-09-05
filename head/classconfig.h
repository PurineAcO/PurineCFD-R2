#pragma once
#include <vector>
#include "config.h"

namespace cc {

struct cell_class;      // 网格
struct face_class;      // 面

struct node_class{
    int number = 0;     // 节点编号
    double x = 0.0, y = 0.0;// 节点坐标
    node_class() = default;
    node_class(int number_,double x_,double y_);// 节点构造器
};

struct face_class{
    int index = 0;      // 面编号
    short type = 0;     // 面类型
    node_class* node[2] = {};      // 面邻接点指针
    vec2 mid = {0.0,0.0};   // 面中点坐标
    vec2 nor = {0.0,0.0};   // 面法向
    int cell_1 = -1, cell_2 = -1;  // 面邻接网格编号
    cell_class* nei[2] = {};       // 面邻接网格指针
    physics phy;                   // 物理量

    face_class() = default;
    face_class(int index_,int p1_,int p2_,int c1_,int c2_,short type_);// 面构造

    // 中心差分重构面上物理量
    void face_physic_mid();  
    // 形成面上所有物理量
    void form_physic();      
};

struct cell_class{
    int index = 0;     // 编号
    int ecnt = 0;      // 面邻接边个数
    int face[4] = {};  // 邻接面编号
    int node[4] = {-1,-1,-1,-1}; // 邻接点编号
    bool fnorm[4] = {}; // 邻接面外法向标记
    face_class* nei[4] = {}; // 邻接面指针
    double vol;         // 体积
    vec2 center;        // 中心坐标

    physics phy;        // 物理量
    double conser[4];   // 守恒量
    double conserformer[4];     // 前期守恒量
    double convect[4];          // 无粘对流项
    dissipation disspiation;    // 耗散项
    double localdt;     // 当地时间步长

    cell_class() = default;
    cell_class(int index_,int f1_,int f2_,int f3_,int f4_,int ecnt_);// 网格构造器

    // 形成全部物理量,用于RK计算开头
    void form_physic();        
    // 形成守恒量
    void form_conservative();  
    // 找到邻接面法向
    void face_normal_out();    
    // 从守恒量解耦物理量,用于RK计算结尾
    void reform();
    // 深度拷贝守恒量,用于时间步
    void copyconver();
};

inline std::vector<node_class> NodeList;    // 节点
inline std::vector<cell_class> CellList;    // 网格
inline std::vector<face_class> FaceList;    // 面
inline std::vector<face_class*> WallFaces;  // 壁面
inline std::vector<face_class*> VILFaces;
inline std::vector<face_class*> POLFaces;
inline std::vector<face_class*> FARFaces;   // 远场

// 安全访问网格
cell_class& gotocell(int number);
// 安全访问邻接边
face_class& gotoface(int number);
// 返回面属性代号
short get_facetype(const char* face_std_name);
// 链接到面
face_class* link_face(int number);
// 链接到网格
cell_class* link_cell(int number);
// 转到邻接网格
cell_class* find_neicell(int index,face_class* face);
// 对于边界边找到内部网格
cell_class* boundary_findcell(face_class* face);

}
