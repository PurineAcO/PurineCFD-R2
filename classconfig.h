#pragma once
#include <utility>
#include <vector>

namespace cc {

struct cell_class;
struct face_class;
struct node_class;

struct physics{

    // 重要变量
    double rho; // 密度
    double u;   // u速度
    double v;   // v速度
    double T;   // 温度

    // 衍生物
    double a;   // 声速
    double mu;  // 空气粘度
    double p;   // 压强

};

struct SA{
    double miubl; // 工作变量miubl
};

struct node_class{
    int number = 0; // 节点编号
    double x = 0.0; // 节点x坐标
    double y = 0.0; // 节点y坐标
    node_class() = default;
    
    // 构造器:记录点的坐标
    node_class(int number_,double x_,double y_);
};

struct face_class{
    int index = 0;                                  // 面编号
    short type = 0;                                 // 面类型
    std::pair<int,int> node = {0,0};           // 节点编号
    std::pair<double,double> mid = {0.0,0.0};  // 面中点坐标
    std::pair<double,double> nor = {0.0,0.0};  // 面法向(未外化)
    int cell_1 = -1; int cell_2 = -1;               // 邻接网格编号(临时使用,常用指针)
    cell_class* nei[2] = {};                        // 邻接网格指针


    face_class() = default;
    // 构造器:形成边,法向量和中点均在此构造完毕.
    face_class(int index_,int p1_,int p2_,int c1_,int c2_,short type_);
};

struct cell_class{
    int index = 0;                                  // 网格索引
    int ecnt = 0;                                   // 网格边数
    int face[4] = {};                               // 网格邻接边编号
    int node[4] = {-1,-1,-1,-1};    // 网格邻接点编号
    face_class* nei[4] = {};                        // 网格邻接边指针
    double vol;                                     // 网格体积
    std::pair<double,double> center ;               // 网格中心
    double sad;                                     // 用于S-A湍流模型的壁面距离
    struct physics phy;                             // 流动变量结构体
    struct SA tur;                                  // 湍流模型字典

    cell_class() = default;
    // 构造器:形成网格
    cell_class(int index_,int f1_,int f2_,int f3_,int f4_,int ecnt_);
};

inline std::vector<node_class> NodeList;     // 节点
inline std::vector<cell_class> CellList;     // 单元
inline std::vector<face_class> FaceList;     // 邻接边
inline std::vector<face_class*> WallFaces;   // 壁面
inline std::vector<face_class*> VILFaces;    // 速度入口
inline std::vector<face_class*> POLFaces;    // 压力出口

}