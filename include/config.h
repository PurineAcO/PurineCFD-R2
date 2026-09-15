#pragma once
#include <string>

namespace cc {

    // 程序基本参数
    inline int cell_num = 0;        // 网格数目
    inline int face_num = 0;        // 面数目
    inline std::string meshpath;    // 网格文件位置
    inline std::string testpath;    // 日志输出路径
    inline std::string fieldpath;   // 流场输出路径
    inline long long max_step = 0;  // 时间步数

    struct ivec2{int x = 0;int y = 0;ivec2() = default;ivec2(int x_,int y_):x(x_),y(y_){}};

    struct vec2{
        double x = 0.0, y = 0.0;
        vec2() = default;
        vec2(double x_,double y_):x(x_),y(y_){}
        vec2& operator+=(const vec2& o){ x += o.x; y += o.y; return *this; }
        vec2& operator-=(const vec2& o){ x -= o.x; y -= o.y; return *this; }
        vec2& operator*=(double s){ x *= s; y *= s; return *this; }
    };

    inline vec2 operator+(vec2 a, const vec2& b){ return a += b; }
    inline vec2 operator-(vec2 a, const vec2& b){ return a -= b; }
    inline vec2 operator*(vec2 a, double s){ return a *= s; }
    inline vec2 operator*(double s, vec2 a){ return a *= s; }
    inline double dot(const vec2& a, const vec2& b){ return a.x*b.x + a.y*b.y; }

    // 物理量矩阵
    struct physics{
        double rho = 0.0, u = 0.0, v = 0.0; // 密度, x/y方向速度
        double T = 0.0, a = 0.0;            // 温度, 声速
        double p = 0.0, e = 0.0;            // 压力, 单位质量总能量
        vec2 ugrad, vgrad, Tgrad;           // 速度与温度梯度
    };

    // 湍流变量矩阵
    struct turbulence{
        double miubl = 0.0;        // SA工作变量
        double miubl_former = 0.0; // 本步RK开始时的miubl
        double miubl_next = 0.0;   // 下一RK阶段的值,算完统一写回
        double sad = 0.0;          // 到最近壁面中点的距离
        vec2 miublgrad;            // miubl梯度
    };

    // 用于JST的人工耗散
    struct dissipation{
        double Y = 0.0;     // 激波捕捉因子
        double L[4] = {};   // 伪Laplace
        double Fd[4] = {};  // JST耗散
    };

    // 支持的边界条件
    inline constexpr short INTER = 0, WALL = 1, FAR = 4;

    // 常数
    inline constexpr double gamma = 1.4;    // 气体绝热常数
    inline constexpr double R = 287.05;     // 气体常数R
    inline constexpr double Cp = 1004.675;  // 定压热容
    inline constexpr double Cv = 717.645;   // 恒容热容
    inline constexpr double Pr = 0.71;      // 普朗特数

}

namespace config {
    bool load(const char* path = "config.json"); // 读取config.json
    inline int dump_step = 1;   // 场输出间隔
    inline int conv_step = 1;   // 残差检查间隔
}
