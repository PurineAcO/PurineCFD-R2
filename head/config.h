#pragma once
#include <cstring>

namespace cc {

    inline double total_time = 0.0;
    inline long long max_step = 0;

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

    struct physics{ double rho, u, v, T, a, p, e; };

    inline constexpr int NEQ = 4;   // 守恒量数

    struct dissipation{
        double Y = 0.0;       // 激波
        double L[NEQ] = {};   // 伪Laplace
        double Fd[NEQ] = {};  // JST
    };

    inline int cell_num = 0, face_num = 0, type_total = 0;
    inline const char* meshpath = "mesh/tunnel.txt";
    inline const char* testpath = "test.txt";
    inline const char* fieldpath = "field";

    inline constexpr short INTER = 0, WALL = 1, POL = 2, VIL = 3, FAR = 4;

    inline struct VIL_condition{double u;double v;double T;double p;} VIL_DEFINE;
    inline struct POL_condition{double p;double T;} POL_DEFINE;
    inline struct VIL_condition FAR_DEFINE;  // 远场自由流

    inline constexpr double gamma = 1.4;
    inline constexpr double R = 287.05;
    inline constexpr double Cp = 1004.675;
    inline constexpr double Cv = 717.645;

}
