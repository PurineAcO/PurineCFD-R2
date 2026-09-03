#pragma once
#include <vector>
#include "config.h"

namespace cc {

struct cell_class;
struct face_class;

struct node_class{
    int number = 0;
    double x = 0.0, y = 0.0;
    node_class() = default;
    node_class(int number_,double x_,double y_);
};

struct face_class{
    int index = 0;
    short type = 0;
    ivec2 node = {0,0};
    vec2 mid = {0.0,0.0};
    vec2 nor = {0.0,0.0};
    int cell_1 = -1, cell_2 = -1;
    cell_class* nei[2] = {};
    physics phy;

    face_class() = default;
    face_class(int index_,int p1_,int p2_,int c1_,int c2_,short type_);

    void face_physic_mid();
    void form_physic();
};

struct cell_class{
    int index = 0;
    int ecnt = 0;
    int face[4] = {};
    int node[4] = {-1,-1,-1,-1};
    bool fnorm[4] = {};
    face_class* nei[4] = {};
    double vol;
    vec2 center;

    physics phy;
    double conser[NEQ];
    double conserformer[NEQ];
    double convect[NEQ];
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

}
