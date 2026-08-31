#pragma once
#include <cstring>
#include <vector>

namespace cc {

    inline double total_time = 0.0;               // 仿真总时间
    inline long long step = 0;                    // 仿真总步数

    struct ivec2{int x = 0;int y = 0;ivec2() = default;ivec2(int x_,int y_):x(x_),y(y_){}};// 二维整数数对
    struct vec2{
        double x = 0.0;
        double y = 0.0;
        vec2() = default;
        vec2(double x_,double y_):x(x_),y(y_){}
        vec2& operator+=(const vec2& o){ x += o.x; y += o.y; return *this; }
        vec2& operator*=(double s)     { x *= s;   y *= s;   return *this; }
    };// 二维double数对

    // 向量加法
    inline vec2 operator+(vec2 a, const vec2& b){ return a += b; }

    // 标量乘法(两侧均可: vec2*s 和 s*vec2)
    inline vec2 operator*(vec2 a, double s){ return a *= s; }
    inline vec2 operator*(double s, vec2 a){ return a *= s; }

    // 点积
    inline double dot(const vec2& a, const vec2& b){ return a.x*b.x + a.y*b.y; }

    struct physics{

    // 重要变量
    double rho; // 密度
    double u;   // u速度
    double v;   // v速度
    double T;   // 温度

    // 不重要物理量
    double a;   // 声速
    double mu;  // 空气粘度
    double p;   // 压强
    double e;   // 比能量

    // 梯度
    vec2 rhograd;   // 密度梯度
    vec2 ugrad;     // u梯度
    vec2 vgrad;     // v梯度
    vec2 Tgrad;     // T梯度

    };// 物理量

    struct SA{
        double miubl;       // 工作变量miubl
        vec2 miublgrad;     // miubl梯度
    };// S-A湍流变量

    struct dissipation{
    double I[3];            // 二阶矩
    double R[2];            // 一阶矩
    vec2 lambda;            // lambda
    double Y;               // 激波探测器
    std::vector<double> L;  // 当前网格的守恒量Laplace算子
    std::vector<double> Fd; // JST阻尼
    };

    // 求解器设置
    inline int cell_num = 0;                          // 网格总数
    inline int face_num = 0;                          // 面总数
    inline int type_total = 0;                        // 边界条件类型总数
    inline const int HALO = 2;                        // HALO网格层数(可能弃用)
    inline const char* meshpath = "mesh/tunnel.txt";  // 网格文件位置
    inline const char* testpath = "test.txt";         // 测试文件位置
    inline const char* turbulence = "SA";             // 湍流模型名称

    using TUR = SA;

    // 网格和边界条件约定
    inline constexpr short INTER = 0;                     // 内部面
    inline constexpr short WALL = 1;                      // 无滑移壁面
    inline constexpr short POL = 2;                       // 压力出口
    inline constexpr short VIL = 3;                       // 速度入口
    inline constexpr short FAR = 4;                       // 压力远场

    // 临时的几何参数
    inline const double H = 2.4;                      // 方腔高度

    // 流动物理条件
    inline struct VIL_condition{double u;double v;            // 速度分量
                                double T;                     // 入口温度,必须是静温
                                double p;                     // 来流静压,亚声速不适用
                                TUR tur;                      // 湍流模型
                                }VIL_DEFINE;   // 速度入口
    inline struct POL_condition{double p;                    // 出口静压
                                double T;                    // 出口温度,必须是总温
                                TUR tur;                      // 湍流模型
                                }POL_DEFINE;   // 压力出口

    // 物理学常数
    inline constexpr double gamma = 1.4;
    inline constexpr double R = 287.05;
    inline constexpr double mu_ref = 1.716e-5;
    inline constexpr double T_ref = 273.15;
    inline constexpr double T_s = 110.4;
    inline constexpr double Cp = 1004.675;
    inline constexpr double Cv = 717.645;

}