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
    node_class(int number_,double x_,double y_);
};

struct face_class{
    int index = 0;      // 面编号
    short type = 0;     // 面类型
    node_class* node[2] = {};         // 面邻接点指针
    vec2 mid = {0.0,0.0};   // 面中点坐标
    vec2 nor = {0.0,0.0};   // 面法向
    int cell_1 = -1, cell_2 = -1;  // 面邻接网格编号
    cell_class* nei[2] = {};       // 面邻接网格指针
    physics phy;                   // 物理量

    face_class() = default;
    face_class(int index_,int p1_,int p2_,int c1_,int c2_,short type_);

    void face_physic_mid();  // 中心差分重构面上物理量
    void form_physic();      // 形成面上所有物理量
};

struct cell_class{
    int index = 0;     // 编号
    int ecnt = 0;      // 面邻接边个数
    int face[4] = {};  // 邻接面编号
    int node[4] = {-1,-1,-1,-1}; // 邻接点指针
    bool fnorm[4] = {}; // 邻接面外法向标记
    face_class* nei[4] = {}; // 邻接面指针
    double vol;         // 体积
    vec2 center;        // 中心坐标

    physics phy;
    double conser[4];
    double conserformer[4];
    double convect[4];
    dissipation disspiation;
    double localdt;

    cell_class() = default;
    cell_class(int index_,int f1_,int f2_,int f3_,int f4_,int ecnt_);

    void form_physic();
    void form_conservative();
    void face_normal_out();
    void reform();
    void copyconver();
};

inline std::vector<node_class> NodeList;
inline std::vector<cell_class> CellList;
inline std::vector<face_class> FaceList;
inline std::vector<face_class*> WallFaces;
inline std::vector<face_class*> VILFaces;
inline std::vector<face_class*> POLFaces;
inline std::vector<face_class*> FARFaces;

cell_class& gotocell(int number);
face_class& gotoface(int number);
short get_facetype(const char* face_std_name);
face_class* link_face(int number);
cell_class* link_cell(int number);
cell_class* find_neicell(int index,face_class* face);
cell_class* boundary_findcell(face_class* face);

}
