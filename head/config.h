#pragma once
#include <cstring>

namespace cc {

    // 程序基本参数
    inline int cell_num = 0;                        // 网格数目
    inline int face_num = 0;                        // 面数目
    inline int type_total = 0;                      // 边界条件数目
    inline const char* meshpath = "mesh/tunnel.txt";// 网格文件位置
    inline const char* testpath = "test.txt";       // 日志输出路径
    inline const char* fieldpath = "field";         // 流场输出路径
    inline double total_time = 0.0;                 // 总时间
    inline long long max_step = 0;                  // 时间步数

    struct ivec2{int x = 0;int y = 0;ivec2() = default;ivec2(int x_,int y_):x(x_),y(y_){}};
    struct vec2{
        double x = 0.0, y = 0.0;
        vec2() = default;
        vec2(double x_,double y_):x(x_),y(y_){}
        vec2& operator+=(const vec2& o){ x += o.x; y += o.y; return *this; }
        vec2& operator*=(double s)     { x *= s;   y *= s;   return *this; }
    };

    inline vec2 operator+(vec2 a, const vec2& b){ return a += b; }
    inline vec2 operator*(vec2 a, double s){ return a *= s; }
    inline vec2 operator*(double s, vec2 a){ return a *= s; }
    inline double dot(const vec2& a, const vec2& b){ return a.x*b.x + a.y*b.y; }

    // 物理量矩阵
    struct physics{ double rho, u, v, T, a, p, e; };

    // 用于JST的人工耗散
    struct dissipation{
        double Y = 0.0;       // 激波捕捉因子
        double L[4] = {};     // 伪Laplace
        double Fd[4] = {};    // JST耗散
    };

    // 支持的边界条件
    inline constexpr short INTER = 0, WALL = 1, POL = 2, VIL = 3, FAR = 4;
    inline struct VIL_condition{double u;double v;double T;double p;} VIL_DEFINE; // 速度入口
    inline struct POL_condition{double p;double T;} POL_DEFINE; // 压力出口
    inline struct VIL_condition FAR_DEFINE;  // 远场自由流

    // 常数
    inline constexpr double gamma = 1.4;    // 气体绝热常数
    inline constexpr double R = 287.05;     // 气体常数R
    inline constexpr double Cp = 1004.675;  // 定压热容
    inline constexpr double Cv = 717.645;   // 恒容热容

}
