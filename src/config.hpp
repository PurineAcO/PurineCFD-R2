#pragma once

#include <cmath>
#include <string>

#define vecfor(cnt) for (int i = 0; i < cnt; i++)

namespace cc {

// 程序基本参数
inline int cell_num = 0;       // 网格数目
inline int face_num = 0;       // 面数目
inline std::string meshpath;   // 网格文件位置
inline std::string testpath;   // 日志输出路径
inline std::string fieldpath;  // 流场输出路径
inline long long max_step = 0; // 时间步数
inline int threads = 0;        // OpenMP线程数, 0表示自动选择
inline bool urans = true;      // 是否是瞬态
inline char scheme = 'R';      // 无粘通量格式: 'R'=Roe, 'J'=JST

struct ivec2 {
    int x = 0;
    int y = 0;
    ivec2() = default;
    ivec2(int x_, int y_) : x(x_), y(y_) {}
};

struct vec2 {
    double x = 0.0, y = 0.0;
    vec2() = default;
    vec2(double x_, double y_) : x(x_), y(y_) {}
    vec2& operator+=(const vec2& o) {
        x += o.x;
        y += o.y;
        return *this;
    }
    vec2& operator-=(const vec2& o) {
        x -= o.x;
        y -= o.y;
        return *this;
    }
    vec2& operator*=(double s) {
        x *= s;
        y *= s;
        return *this;
    }
    void clear();
    double norm();
};

inline vec2 operator+(vec2 a, const vec2& b) {
    return a += b;
}
inline vec2 operator-(vec2 a, const vec2& b) {
    return a -= b;
}
inline vec2 operator*(vec2 a, double s) {
    return a *= s;
}
inline vec2 operator*(double s, vec2 a) {
    return a *= s;
}
inline double dot(const vec2& a, const vec2& b) {
    return a.x * b.x + a.y * b.y;
}
inline void vec2::clear() {
    x = 0.0;
    y = 0.0;
}
inline double vec2::norm() {
    return std::sqrt(x * x + y * y);
}

struct vec4 {
    double c; // continuous
    double x; // velocity-x
    double y; // velocity-y
    double e; // energy
    vec4() = default;
    vec4(double c_, double x_, double y_, double e_) : c(c_), x(x_), y(y_), e(e_) {}
    vec4& operator+=(const vec4& o) {
        c += o.c;
        x += o.x;
        y += o.y;
        e += o.e;
        return *this;
    }
    vec4& operator-=(const vec4& o) {
        c -= o.c;
        x -= o.x;
        y -= o.y;
        e -= o.e;
        return *this;
    }
    vec4& operator*=(double k) {
        c *= k;
        x *= k;
        y *= k;
        e *= k;
        return *this;
    }
    void clear();
};

inline vec4 operator+(vec4 a, const vec4& b) {
    return a += b;
}
inline vec4 operator-(vec4 a, const vec4& b) {
    return a -= b;
}
inline vec4 operator*(vec4 a, double s) {
    return a *= s;
}
inline vec4 operator*(double s, vec4 a) {
    return a *= s;
}
inline void vec4::clear() {
    c = 0.0;
    x = 0.0;
    y = 0.0;
    e = 0.0;
}

struct vecp {
    double rho, u, v, T;
    vecp() = default;
    vecp(double rho_, double u_, double v_, double T_) : rho(rho_), u(u_), v(v_), T(T_) {}
    vecp& operator+=(const vecp& o) {
        rho += o.rho;
        u += o.u;
        v += o.v;
        T += o.T;
        return *this;
    }
    vecp& operator-=(const vecp& o) {
        rho -= o.rho;
        u -= o.u;
        v -= o.v;
        T -= o.T;
        return *this;
    }
    vecp& operator*=(double k) {
        rho *= k;
        u *= k;
        v *= k;
        T *= k;
        return *this;
    }
    void clear();
};

inline vecp operator+(vecp a, const vecp& b) {
    return a += b;
}
inline vecp operator-(vecp a, const vecp& b) {
    return a -= b;
}
inline vecp operator*(vecp a, double s) {
    return a *= s;
}
inline vecp operator*(double s, vecp a) {
    return a *= s;
}
inline void vecp::clear() {
    rho = 0.0;
    u = 0.0;
    v = 0.0;
    T = 0.0;
}

struct vecgrad {
    vec2 rhograd, ugrad, vgrad, Tgrad;
    vecgrad() = default;
    vecgrad(vec2 rho_, vec2 u_, vec2 v_, vec2 t_)
        : rhograd(rho_), ugrad(u_), vgrad(v_), Tgrad(t_) {}
    void clear();
};
inline vecp operator*(const vecgrad& a, const vec2& r) {
    return vecp(dot(a.rhograd, r), dot(a.ugrad, r), dot(a.vgrad, r), dot(a.Tgrad, r));
}
inline vecgrad operator*(const vecp& phyf, const vec2& fnor) {
    return vecgrad(phyf.rho * fnor, phyf.u * fnor, phyf.v * fnor, phyf.T * fnor);
}
inline void vecgrad::clear() {
    rhograd.clear();
    ugrad.clear();
    vgrad.clear();
    Tgrad.clear();
}
inline vecgrad operator*(const double a, const vecgrad& grad) {
    return vecgrad(a * grad.rhograd, a * grad.ugrad, a * grad.vgrad, a * grad.Tgrad);
}
inline vecgrad operator+(const vecgrad& a, const vecgrad& b) {
    return vecgrad(a.rhograd + b.rhograd, a.ugrad + b.ugrad, a.vgrad + b.vgrad, a.Tgrad + b.Tgrad);
}

struct otphy {
    double a, p, e, un, mu;
    // 从原始物理量生成引申物理量
    otphy() = default;
    otphy(double a_, double p_, double e_) : a(a_), p(p_), e(e_) {}
    void form_otphy(cc::vecp phy);
};

// 湍流变量矩阵
struct turbulence {
    double miubl = 0.0;        // SA工作变量
    double miubl_former = 0.0; // 本步RK开始时的miubl
    double miubl_next = 0.0;   // 下一RK阶段的值,算完统一写回
    double sad = 0.0;          // 到最近壁面中点的距离
    vec2 miublgrad;            // miubl梯度
    double mueff = 0.0;        // 有效粘性系数
};

// 用于JST的人工耗散
struct dissipation {
    double Y = 0.0;    // 激波捕捉因子
    double L[4] = {};  // 伪Laplace
    double Fd[4] = {}; // JST耗散
};

// 支持的边界条件
inline constexpr short INTER = 0, WALL = 1, FAR = 4;

// 常数
inline constexpr double gamma = 1.4;   // 气体绝热常数
inline constexpr double R = 287.05;    // 气体常数R
inline constexpr double Cp = 1004.675; // 定压热容
inline constexpr double Cv = 717.645;  // 恒容热容
inline constexpr double Pr = 0.71;     // 普朗特数

} // namespace cc

// 伪时间步参数
namespace fatime {
inline double CFL = 1.0;
}

// 求解设置
namespace config {
inline bool load(const char* path = "config.json"); // 读取config.json, 实现在 io.hpp
inline int dump_step = 1;                           // 场输出间隔
inline int conv_step = 1;                           // 残差检查间隔
} // namespace config

// URANS双时间步参数, 对应config.json中的solver.urans
namespace urans {
inline double dt = 0.0;          // 物理时间步长, s
inline int inner = 20;           // 每个物理时间步的最大内迭代次数
inline double inner_tol = 1e-3;  // 内迭代残差相对首次内迭代的下降目标
inline int sweeps = 4;           // 每次内迭代的对称Gauss-Seidel扫描次数
inline int steady_iters = 0;     // URANS之前的定常隐式迭代次数, 用来给初场
inline double steady_cfl = 50.0; // 定常迭代的CFL终值(从2线性增大到该值)
inline int wall_interval = 0;    // 壁面Cp输出间隔, 0表示不输出
inline double seed = 0.0;        // 初场反对称涡扰动幅值(相对U∞), 对称问题起振用
} // namespace urans

// 结构化网格参数
namespace structer {
inline bool ifstructer = false; // 结构化网格令牌
inline std::string adjacency;   // 结构化邻接表路径
inline int S_MAX = 0;           // 环向单元数, 由邻接表表头给出
inline int N_MAX = 0;           // 径向单元数, 由邻接表表头给出
inline constexpr int HALO = 3;  // HALO网格层数
} // namespace structer
