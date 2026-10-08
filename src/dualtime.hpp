#pragma once

#include "HALO.hpp"
#include "SA.hpp"
#include "boundary.hpp"
#include "classconfig.hpp"
#include "config.hpp"
#include "convect.hpp"
#include "grad.hpp"
#include "interpolate.hpp"
#include <algorithm>
#include <cmath>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

/*
双时间步URANS, 做法与Fluent密度基隐式求解器相同:
  外迭代(物理时间)用一阶后向欧拉:  V(Q^{n+1}-Q^n)/Δt + R(Q^{n+1}) = 0
  内迭代(伪时间)隐式推进, 每步解线性方程
      [V/Δτ + V/Δt + ∂R/∂Q] ΔQ = -[R(Q^k) + V(Q^k-Q^n)/Δt]
  ∂R/∂Q 用一阶Rusanov通量的近似Jacobian:
      对角  0.5*Σλ_f (标量)
      邻居j 0.5*(F(Q_j+ΔQ_j)-F(Q_j))·n - 0.5*λ_f*ΔQ_j
  λ_f 为面上对流+黏性谱半径. 线性方程用对称Gauss-Seidel扫描求解,
  线程按单元编号分块, 块内用最新值, 块间用上一遍的值, 结果与线程调度无关.
  SA方程在流动方程之后单独求解(与Fluent一致), 同样BDF1+隐式伪时间,
  对流项一阶迎风, 扩散项两点格式, 破坏项线性化后放进对角.
  dt<=0时去掉物理时间项, 即隐式定常迭代.
*/

namespace dual {
// 计算整场残差 R = convect - visflux (单元净流出通量)
void eval_residual();
// 流动方程一次隐式伪时间步, 返回连续方程残差的L2
double flow_step(double dt, double cfl);
// SA方程一次隐式伪时间步, 需紧跟flow_step调用(复用其面通量与梯度)
void sa_step(double dt, double cfl);
} // namespace dual

namespace dual {

inline std::vector<cc::vec4> rhs, dQ, dQold;
inline std::vector<double> lam, diag, nrhs, dnu, dnuold;

inline void eval_residual() {
#pragma omp parallel for schedule(static)
    allcell icell(i).prim();
    noslip_wall_boundary();
    far_field_boundary();
    if (structer::ifstructer) update_ghost_field();
#pragma omp parallel for schedule(static)
    allfac {
        interpolate_mid(&iface(i));
        muscl(&iface(i));
    }
#pragma omp parallel for schedule(static)
    allcell least_square_cell_based(icell(i));
#pragma omp parallel for schedule(static)
    allcell grad_onface(icell(i));
#pragma omp parallel for schedule(static)
    allfac convect_ROE(iface(i));
#pragma omp parallel for schedule(static)
    allcell assemble_flux(icell(i), 'R');
#pragma omp parallel for schedule(static)
    allfac SA::diffusion_SA(iface(i));
#pragma omp parallel for schedule(static)
    allcell SA::assemble_visflux(icell(i));
}

// 守恒量Q沿n(带面长)的无粘通量
static inline cc::vec4 euler_flux(const cc::vec4& Q, double nx, double ny) {
    const double u = Q.x / Q.c, v = Q.y / Q.c;
    const double p = (cc::gamma - 1.0) * (Q.e - 0.5 * Q.c * (u * u + v * v));
    const double un = u * nx + v * ny;
    return cc::vec4(Q.c * un, Q.x * un + p * nx, Q.y * un + p * ny, (Q.e + p) * un);
}

// 面谱半径(带面长): 对流(|un|+a)S + 黏性 μeff/ρ*S²/V
static inline double face_lambda(const cc::face_class* f, const cc::cell_class& cell) {
    const double un = std::abs(f->phy.u * f->nor.x + f->phy.v * f->nor.y) / f->len;
    const double Cn = std::max(4.0 / 3.0, cc::gamma / cc::Pr);
    return (un + f->otphy.a) * f->len + Cn * f->tur.mueff / f->phy.rho * f->len * f->len / cell.vol;
}

// 内部面另一侧的单元, 边界面返回nullptr
static inline cc::cell_class* other_cell(const cc::face_class* f, const cc::cell_class& cell) {
    if (f->type != cc::INTER) return nullptr;
    return f->nei[0] == &cell ? f->nei[1] : f->nei[0];
}

// 分块对称GS: 每个线程先正扫再反扫自己的单元区间[lo,hi)
template <class Relax> static void block_sgs(Relax&& relax) {
#pragma omp parallel
    {
#ifdef _OPENMP
        const int nt = omp_get_num_threads(), tid = omp_get_thread_num();
#else
        const int nt = 1, tid = 0;
#endif
        const int lo = (int)((long long)cc::cell_num * tid / nt);
        const int hi = (int)((long long)cc::cell_num * (tid + 1) / nt);
        for (int c = lo; c < hi; c++)
            relax(c, lo, hi);
        for (int c = hi - 1; c >= lo; c--)
            relax(c, lo, hi);
    }
}

inline double flow_step(double dt, double cfl) {
    const int N = cc::cell_num;
    if ((int)rhs.size() != N) {
        rhs.resize(N);
        dQ.resize(N);
        dQold.resize(N);
        lam.resize(N);
        diag.resize(N);
        nrhs.resize(N);
        dnu.resize(N);
        dnuold.resize(N);
    }
    eval_residual();
    const double inv_dt = dt > 0.0 ? 1.0 / dt : 0.0;
    double res2 = 0.0;
#pragma omp parallel for schedule(static) reduction(+ : res2)
    for (int c = 0; c < N; c++) {
        cc::cell_class& cell = icell(c);
        double sum = 0.0;
        allface(cell) sum += face_lambda(cell.faces[i], cell);
        lam[c] = sum;
        // V/Δτ = Σλ/CFL
        diag[c] = cell.vol * inv_dt + sum / cfl + 0.5 * sum;
        rhs[c] = cell.convect - cell.visflux;
        if (dt > 0.0) rhs[c] += (cell.vol * inv_dt) * (cell.conser - cell.lastconser);
        dQ[c] = cc::vec4(0, 0, 0, 0);
        res2 += rhs[c].c * rhs[c].c;
    }
    auto relax = [](int c, int lo, int hi) {
        const cc::cell_class& cell = icell(c);
        cc::vec4 acc = rhs[c];
        allface(cell) {
            const cc::face_class* f = cell.faces[i];
            const cc::cell_class* nb = other_cell(f, cell);
            if (!nb) continue;
            const int j = nb->index - 1;
            const cc::vec4& dj = (j >= lo && j < hi) ? dQ[j] : dQold[j];
            // 本单元的外法向
            const double s = 2 * cell.fnorm[i] - 1;
            const double nx = s * f->nor.x, ny = s * f->nor.y;
            acc += 0.5 * (euler_flux(nb->conser + dj, nx, ny) - euler_flux(nb->conser, nx, ny)) -
                   (0.5 * face_lambda(f, cell)) * dj;
        }
        dQ[c] = acc * (-1.0 / diag[c]);
    };
    for (int sw = 0; sw < urans::sweeps; sw++) {
#pragma omp parallel for schedule(static)
        for (int c = 0; c < N; c++)
            dQold[c] = dQ[c];
        block_sgs(relax);
    }
// 更新守恒量; 若密度或内能下降超过一半则把该单元的增量减半, 保证正性
#pragma omp parallel for schedule(static)
    for (int c = 0; c < N; c++) {
        const cc::vec4 Q = icell(c).conser;
        const double ein = Q.e - 0.5 * (Q.x * Q.x + Q.y * Q.y) / Q.c;
        double w = 1.0;
        for (int it = 0; it < 6; it++) {
            const cc::vec4 Qn = Q + w * dQ[c];
            if (Qn.c > 0.5 * Q.c && Qn.e - 0.5 * (Qn.x * Qn.x + Qn.y * Qn.y) / Qn.c > 0.5 * ein)
                break;
            w *= 0.5;
        }
        icell(c).conser = Q + w * dQ[c];
    }
    return std::sqrt(res2 / N);
}

// SA破坏项对ν̃的导数, 近似为 2ρ*Cw1*fw*ν̃/d² (fw视为常数)
static inline double sa_destruction_jac(const cc::cell_class& cell) {
    const double nu = std::max(cell.tur.miubl, 0.0);
    const double chi = cell.phy.rho * nu / cell.otphy.mu;
    const double fv1 = chi * chi * chi / (chi * chi * chi + SA::Cv1 * SA::Cv1 * SA::Cv1);
    const double fv2 = 1 - chi / (1 + chi * fv1);
    const double d2 = cell.tur.sad * cell.tur.sad;
    const double vort = std::abs(cell.phgrad.vgrad.x - cell.phgrad.ugrad.y);
    const double sn = nu / (d2 * SA::kappa * SA::kappa);
    const double St = std::max(vort + fv2 * sn, std::max(0.3 * vort, 1e-20));
    const double r = std::min(sn / St, SA::rmax);
    const double g = r + SA::Cw2 * (std::pow(r, 6) - r);
    const double c6 = std::pow(SA::Cw3, 6);
    const double fw = g * std::pow((1 + c6) / (std::pow(g, 6) + c6), 1.0 / 6);
    return 2.0 * cell.phy.rho * SA::Cw1 * fw * nu / d2;
}

// 面的SA扩散系数 (μ+ρν̃)/σ
static inline double sa_coef(const cc::face_class* f) {
    return SA::inv_sigma * (f->otphy.mu + f->phy.rho * std::max(f->tur.miubl, 0.0));
}

// 单元中心到面另一侧(邻居中心或面中点)沿面法向的距离
static inline double normal_dist(const cc::face_class* f, const cc::cell_class& cell,
                                 const cc::cell_class* nb) {
    const cc::vec2 d = nb ? nb->center - cell.center : f->mid - cell.center;
    return std::abs(cc::dot(d, f->nor)) / f->len;
}

inline void sa_step(double dt, double cfl) {
    const int N = cc::cell_num;
    const double inv_dt = dt > 0.0 ? 1.0 / dt : 0.0;
// 守恒形式 ∂(ρν̃)/∂t + Σ(ṁν̃) - Σ(coef ∂ν̃/∂n S) - V*S_SA = 0, 未知量为ν̃
#pragma omp parallel for schedule(static)
    for (int c = 0; c < N; c++) {
        const cc::cell_class& cell = icell(c);
        const double nu = cell.tur.miubl;
        double res = 0.0, dg = 0.0;
        allface(cell) {
            const cc::face_class* f = cell.faces[i];
            const cc::cell_class* nb = other_cell(f, cell);
            // Roe质量通量沿nei[0]->nei[1], 换成本单元外流
            const double mdot = (f->nei[0] == &cell ? 1.0 : -1.0) * f->convect.c;
            const double nu_o = nb ? nb->tur.miubl : f->tur.miubl;
            res += mdot * (mdot > 0.0 ? nu : nu_o);
            const double coef = sa_coef(f);
            double dist = 0.0, dnudn = 0.0;
            if (nb) {
                // 插值面梯度沿中心连线方向用两点差分修正, 避免奇偶失耦
                const cc::vec2 d = nb->center - cell.center;
                const double len = std::sqrt(cc::dot(d, d));
                const cc::vec2 e = d * (1.0 / len);
                const cc::vec2 g = f->tur.miublgrad;
                const double corr = (nu_o - nu) / len - cc::dot(g, e);
                dnudn = cc::dot(g + corr * e, (2 * cell.fnorm[i] - 1) * f->nor);
                dist = normal_dist(f, cell, nb);
            } else {
                // 壁面ν̃=0、距离取壁面距离; 远场取面值
                dist = (f->type == cc::WALL) ? cell.tur.sad : normal_dist(f, cell, nullptr);
                dnudn = (nu_o - nu) / dist * f->len;
            }
            res -= coef * dnudn;
            dg += std::max(mdot, 0.0) + coef * f->len / dist;
        }
        res -= cell.vol * SA::source_SA(cell);
        if (dt > 0.0)
            res += cell.vol * inv_dt *
                   (cell.conser.c * nu - cell.lastconser.c * cell.tur.miubl_former);
        dg += cell.conser.c * (cell.vol * inv_dt + lam[c] / cfl) +
              cell.vol * sa_destruction_jac(cell);
        nrhs[c] = res;
        diag[c] = dg;
        dnu[c] = 0.0;
    }
    auto relax = [](int c, int lo, int hi) {
        const cc::cell_class& cell = icell(c);
        double acc = nrhs[c];
        allface(cell) {
            const cc::face_class* f = cell.faces[i];
            const cc::cell_class* nb = other_cell(f, cell);
            if (!nb) continue;
            const int j = nb->index - 1;
            const double dj = (j >= lo && j < hi) ? dnu[j] : dnuold[j];
            const double mdot = (f->nei[0] == &cell ? 1.0 : -1.0) * f->convect.c;
            acc += (std::min(mdot, 0.0) - sa_coef(f) * f->len / normal_dist(f, cell, nb)) * dj;
        }
        dnu[c] = -acc / diag[c];
    };
    for (int sw = 0; sw < urans::sweeps; sw++) {
#pragma omp parallel for schedule(static)
        for (int c = 0; c < N; c++)
            dnuold[c] = dnu[c];
        block_sgs(relax);
    }
#pragma omp parallel for schedule(static)
    for (int c = 0; c < N; c++) {
        double dn = dnu[c];
        // 只丢非有限值; 不能做幅值限幅, 否则 ν̃ 从 0 起长时会被钉死
        if (!std::isfinite(dn)) { dn = 0.0; }
        icell(c).tur.miubl = std::max(0.0, icell(c).tur.miubl + dn);
    }
}

} // namespace dual
