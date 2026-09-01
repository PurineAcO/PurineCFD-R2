#pragma once
#include <vector>
#include "config.h"

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
    int index = 0;                                 // 面编号
    short type = 0;                                // 面类型
    ivec2 node = {0,0};                     // 节点编号
    vec2 mid = {0.0,0.0};                   // 面中点坐标
    vec2 nor = {0.0,0.0};                   // 面法向(未外化)
    int cell_1 = -1; int cell_2 = -1;              // 邻接网格编号(临时使用,常用指针)
    cell_class* nei[2] = {};                       // 邻接网格指针
    physics phy;
    TUR tur;

    
    // 构造器:形成边,法向量和中点均在此构造完毕.
    face_class() = default;
    face_class(int index_,int p1_,int p2_,int c1_,int c2_,short type_);

    // 中心差分计算面上物理量
    void face_physic_mid();
    // 形成全部物理量
    void form_physic();
};

struct cell_class{
    
    // 基本
    int index = 0;                                  // 网格索引
    int ecnt = 0;                                   // 网格边数
    int face[4] = {};                               // 网格邻接边编号
    int node[4] = {-1,-1,-1,-1};    // 网格邻接点编号
    bool fnorm[4] = {};                             // 网格邻接边法向量标记
    face_class* nei[4] = {};                        // 网格邻接边指针
    double vol;                                     // 网格体积
    vec2 center;                                    // 网格中心

    // 物理
    double sad;                                     // 用于S-A湍流模型的壁面距离
    physics phy;                                    // 流动变量结构体
    TUR tur;                                        // 湍流模型字典
    std::vector<double> conser;                     // RK上的守恒量
    std::vector<double> conserformer;               // 时间步上的守恒量
    std::vector<double> convect;                    // 对流项
    std::vector<double> diffusion;                  // 扩散项
    std::vector<double> source;                     // 源项
    dissipation disspiation;                        // 耗散
    double localdt;                                 // 当地时间步长

    // 构造器:形成网格
    cell_class() = default;
    cell_class(int index_,int f1_,int f2_,int f3_,int f4_,int ecnt_);

    // 形成全部物理量
    void form_physic();
    // 形成守恒量
    void form_conservative();
    // 判断邻接面的法向
    void face_normal_out();
    // 还原物理量
    void reform();
    // 从RK守恒量转到时间步上守恒量
    void copyconver();
};

inline std::vector<node_class> NodeList;     // 节点
inline std::vector<cell_class> CellList;     // 单元
inline std::vector<face_class> FaceList;     // 邻接边
inline std::vector<face_class*> WallFaces;   // 壁面
inline std::vector<face_class*> VILFaces;    // 速度入口
inline std::vector<face_class*> POLFaces;    // 压力出口

// 安全的网格访问
cell_class& gotocell(int number);

// 安全的面访问
face_class& gotoface(int number);

// 判断面的类型
short get_facetype(const char* face_std_name);

// 链接一个面
face_class* link_face(int number);

// 链接一个网格
cell_class* link_cell(int number);

// 找到邻接网格
cell_class* find_neicell(int index,face_class* face);

}