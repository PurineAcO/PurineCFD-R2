#pragma once

#include "SA.hpp"
#include "classconfig.hpp"
#include "config.hpp"
#include <algorithm>
#include <cmath>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace dual {
// SA方程一次隐式伪时间步, 需紧跟 impl::flow_step 调用(复用其面通量与隐式谱半径)
void sa_step(double dt, double cfl);
} // namespace dual

namespace dual {

inline std::vector<double> diag, nrhs, dnu, dnuold;

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
    // 本步不再依赖 dual::flow_step 预分配, 自己兜底
    if ((int)nrhs.size() != N) {
        nrhs.resize(N);
        diag.resize(N);
        dnu.resize(N);
        dnuold.resize(N);
    }
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
        dg += cell.conser.c * (cell.vol * inv_dt + cell.lambda / cfl) +
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
