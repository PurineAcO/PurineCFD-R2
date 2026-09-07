#pragma once

namespace cfd {

struct Vector2 {
  double x = 0.0, y = 0.0;
  Vector2() = default;
  Vector2(double x_, double y_) : x(x_), y(y_) {}
  Vector2& operator+=(const Vector2& o) {
    x += o.x;
    y += o.y;
    return *this;
  }
  Vector2& operator*=(double s) {
    x *= s;
    y *= s;
    return *this;
  }
};

inline Vector2 operator+(Vector2 a, const Vector2& b) {
  return a += b;
}
inline Vector2 operator*(Vector2 a, double s) {
  return a *= s;
}
inline Vector2 operator*(double s, Vector2 a) {
  return a *= s;
}
inline double dot(const Vector2& a, const Vector2& b) {
  return a.x * b.x + a.y * b.y;
}

// 原始变量采用 SI 单位；e 是单位质量总能量，不是内能。
struct FlowState {
  double rho = 0.0;        // 密度，kg/m³
  double u = 0.0, v = 0.0; // x、y 方向速度，m/s
  double T = 0.0;          // 温度，K
  double a = 0.0;          // 声速，m/s
  double p = 0.0;          // 压力，Pa
  double e = 0.0;          // E = Cv*T + (u²+v²)/2，J/kg
  Vector2 ugrad, vgrad, Tgrad;
};

// Spalart–Allmaras 工作变量 ν̃；它不直接等于湍流运动黏度 νt。
struct TurbulenceState {
  double nu_tilde = 0.0;                      // m²/s；νt = ν̃*fv1
  double nu_tilde_previous = 0.0;             // 本步 RK 开始时的值
  double nu_tilde_next = 0.0;                 // 下一 RK 阶段的值，统一同步后写回
  double inverse_wall_distance_squared = 0.0; // 固定几何量 1/d²
  Vector2 nu_tilde_gradient;
};

// 用于JST的人工耗散
struct JstDissipation {
  double pressure_sensor = 0.0; // 激波捕捉因子
  double laplacian[4] = {};     // 伪Laplace
  double flux[4] = {};          // JST耗散
};

struct Freestream {
  double u = 0.0, v = 0.0, T = 0.0, p = 0.0;
};
inline Freestream freestream;

// 理想空气常数；R、Cp、Cv 的单位为 J/(kg·K)。
inline constexpr double gamma = 1.4;   // 气体绝热常数
inline constexpr double R = 287.05;    // 气体常数R
inline constexpr double Cp = 1004.675; // 定压热容
inline constexpr double Cv = 717.645;  // 恒容热容
inline constexpr double Pr = 0.71;     // 普朗特数

} // namespace cfd
