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

// 流动变量采用 SI 单位。
struct FlowState {
  double rho = 0.0;        // 密度，kg/m³
  double u = 0.0, v = 0.0; // x、y 方向速度，m/s
  double T = 0.0;          // 温度，K
  double a = 0.0;          // 声速，m/s
  double p = 0.0;          // 压力，Pa
  double e = 0.0;          // 单位质量总能量 E = Cv*T + (u²+v²)/2，J/kg
  Vector2 ugrad, vgrad, Tgrad;
};

// Spalart-Allmaras 模型用工作变量 ν̃ 计算湍流运动黏度 νt = ν̃*fv1。
struct TurbulenceState {
  double nu_tilde = 0.0;                      // SA 工作变量，m²/s
  double nu_tilde_previous = 0.0;             // 本步 RK 开始时的值
  double nu_tilde_next = 0.0;                 // 下一 RK 阶段的值，统一同步后写回
  double inverse_wall_distance_squared = 0.0; // 固定几何量 1/d²
  Vector2 nu_tilde_gradient;
};

// JST 人工耗散所需的压力传感器、守恒量差分和耗散通量。
struct JstDissipation {
  double pressure_sensor = 0.0; // 压力传感器 Σⱼ|pⱼ-pᵢ| / Σⱼ(pⱼ+pᵢ)
  double laplacian[4] = {};     // 守恒量差分 Lᵢ = Σⱼ(Qⱼ-Qᵢ)
  double flux[4] = {};          // 单元各面的 JST 耗散通量之和
};

struct Freestream {
  double u = 0.0, v = 0.0, T = 0.0, p = 0.0;
};
inline Freestream freestream;

// 理想空气常数；R、Cp、Cv 的单位为 J/(kg·K)。
inline constexpr double gamma = 1.4;   // 比热比 Cp/Cv
inline constexpr double R = 287.05;    // 空气比气体常数
inline constexpr double Cp = 1004.675; // 定压比热容
inline constexpr double Cv = 717.645;  // 定容比热容
inline constexpr double Pr = 0.71;     // 普朗特数 Pr

} // namespace cfd
