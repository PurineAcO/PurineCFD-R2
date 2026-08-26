#pragma once
#include <vector>

namespace cc {

struct cell_class;
struct face_class;
struct node_class;

struct node_class{
    int number = 0; // 节点编号
    double x = 0.0; // 节点x坐标
    double y = 0.0; // 节点y坐标
    node_class() = default;
    
    // 构造器:记录点的坐标
    node_class(int number_,double x_,double y_);
};

struct face_class{
    int index = 0; // 面编号
    std::pair<double,double> mid = {0.0,0.0};
    std::pair<double,double> nor = {0.0,0.0};
    int cell_1 = -1; int cell_2 = -1;
    face_class() = default;

    // 构造器:形成边,法向量和中点均在此构造完毕.
    face_class(int index_,int p1_,int p2_,int c1_,int c2_);
};

struct cell_class{
    int index = 0;int ecnt = 0;int face[4] = {};
    cell_class() = default;

    // 构造器:形成网格
    cell_class(int index_,int f1_,int f2_,int f3_,int f4_,int ecnt_);
};

inline std::vector<node_class> NodeList;
inline std::vector<cell_class> CellList;
inline std::vector<face_class> FaceList;

}