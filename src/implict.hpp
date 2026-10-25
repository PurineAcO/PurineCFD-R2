#pragma once

#include "HALO.hpp"
#include "SA.hpp"
#include "boundary.hpp"
#include "classconfig.hpp"
#include "config.hpp"
#include "convect.hpp"
#include "grad.hpp"
#include "interpolate.hpp"
#include <cmath>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace impl {

// 松弛单元数
inline int n_clip = 0;

// 矩阵量
inline std::vector<cc::vec4> dQ;   // 守恒量变化量
inline std::vector<cc::vec4> dQold;// 冻结守恒量变化量
inline std::vector<cc::vec4> rhs;  // 右侧残差

// 形成右侧残差
void right_residual();
// 分配求解器向量
void malloc_vector();
// 清理矩阵
void clear_mat();
// 形成无粘通量Jacobi
void form_euler_jacobi(cc::cell_class& cell);
// 形成粘性通量Jacobi
void form_viscous_jacobi(cc::face_class& face);
// 壁面特别处理
void wall_jacobi(cc::face_class *face,cc::cell_class* cell,double s);
// 远场特别处理
void far_jacobi(cc::face_class *face,cc::cell_class* cell,double s);
// 伪时间谱半径
double face_lambda(const cc::face_class* f, const cc::cell_class& cell);
// 求出对角线
void form_diag(cc::cell_class& cell,double inv_dt,double cfl);
// 复制解向量到缓存
void Qcopy();
// 求解核心
void solvekernel(int c,int lo,int hi);
// 循环器GS
template <class KERNEL> void block_sgs(KERNEL&& solve_kernel);
// 松弛
int relax();
// 进行一次双隐式伪时间迭代
double flow_step(double dt, double cfl);
// 进行一次湍流模型伪时间迭代
void turbulence_step(double dt, double cfl);

}// namespace impl

inline void impl::right_residual(){
#pragma omp parallel for schedule(static)
    allcell icell(i).prim();
    noslip_wall_boundary();
    far_field_boundary();
    if(structer::ifstructer) update_ghost_field();
#pragma omp parallel for schedule(static)
    allfac{interpolate_mid(&iface(i));muscl(&iface(i));}
#pragma omp parallel for schedule(static)
    allcell least_square_cell_based(icell(i));
#pragma omp parallel for schedule(static)
    allcell grad_onface(icell(i));
#pragma omp parallel for schedule(static)
    allfac convect_ROE(iface(i));
#pragma omp parallel for schedule(static)
    allcell assemble_flux(icell(i),'R');
#pragma omp parallel for schedule(static)
    allfac SA::diffusion_SA(iface(i));
#pragma omp parallel for schedule(static)
    allcell SA::assemble_visflux(icell(i));
}

inline void impl::malloc_vector(){
    if((int) dQ.size() == cc::cell_num) return;
    dQ.resize(cc::cell_num);
    dQold.resize(cc::cell_num);
    rhs.resize(cc::cell_num);
}

inline void impl::clear_mat(){
#pragma omp parallel for schedule(static)
    allcell icell(i).diag.clear();
#pragma omp parallel for schedule(static)
    allcell icell(i).lambda = 0.0;
    // SA::source_jacobi_clear();
}

inline void impl::form_euler_jacobi(cc::cell_class& cell){
    const std::pair<cc::mat5,cc::mat5> m = form_pfgpw_convect(cell);
    cell.conj = std::pair<cc::mat4, cc::mat4>(cc::slice4(m.first),cc::slice4(m.second));
}

inline void impl::form_viscous_jacobi(cc::face_class & face){
    const std::pair<cc::mat5,cc::mat5> dm = SA::form_pfpgq_diffusion(face);
    face.imats.Bf = cc::slice4(face.toface_jacobi(std::pair<cc::mat5,cc::mat5>(
            SA::form_pfpw_from_pfpgq(dm.first,face),
            SA::form_pfpw_from_pfpgq(dm.second,face))));
    if(face.type != cc::INTER){
        face.imats.Ab.clear();
        face.imats.Tf.clear();
        return;
    }
    face.imats.Ab = cc::slice4(roe_abs_n(face));
    face.imats.Tf = cc::slice4(SA::diffusion_grad_jac(face,normal_dist(&face)));
}

inline void impl::wall_jacobi(cc::face_class *face,cc::cell_class* cell,double s){
    // 发生在组装期
    const double lamc = (std::abs(face->phy.u*face->nor.x + face->phy.v*face->nor.y)/face->len+face->otphy.a) * face->len;
    cell->diag += 0.5*lamc*cc::eyes;
    cell->diag += cc::slice4(SA::diffusion_grad_jac(*face, cell->tur.sad));
    cell->diag += (-s)*face->imats.Bf;
}

inline void impl::far_jacobi(cc::face_class *face, cc::cell_class *cell, double s){
    // 发生在组装期
    cell->diag += (-s)*face->imats.Bf;
}

inline double impl::face_lambda(const cc::face_class* f, const cc::cell_class& cell) {
    const double un = std::abs(f->phy.u * f->nor.x + f->phy.v * f->nor.y) / f->len;
    const double Cn = std::max(4.0 / 3.0, cc::gamma / cc::Pr);
    return (un + f->otphy.a) * f->len + Cn * f->tur.mueff / f->phy.rho * f->len * f->len / cell.vol;
}

inline void impl::form_diag(cc::cell_class &cell,double inv_dt,double cfl){
    // 格心源对角
    allface(cell) cell.lambda += face_lambda(cell.faces[i], cell);
    cell.diag += (cell.vol*inv_dt + cell.lambda/cfl) * cc::eyes;
    cc::mat4 far;far.clear();
    // 面上的对角部分
    allface(cell){
        double s = 2*cell.fnorm[i]-1;
        cc::face_class* f = cell.faces[i];
        if(f->type == cc::INTER) cell.diag += 0.5*s*f->toface_jacobi(cell.conj) + 0.5*f->imats.Ab +f->imats.Tf + (-0.5*s)*f->imats.Bf;
        else if(f->type == cc::WALL) wall_jacobi(f, &cell, s);
        else if(f->type == cc::FAR){ far_jacobi(f,&cell, s); far += far_convect_jacobi(f); }
    }
    cell.diag += far;
}

inline void impl::Qcopy(){
#pragma omp parallel for schedule(static)
    for(int i=0;i<(int) dQ.size();i++) dQold[i]=dQ[i];
}

inline void impl::solvekernel(int c, int lo, int hi){
    const cc::cell_class& cell = icell(c);
    cc::vec4 acc = rhs[c]; 
    allface(cell){
        const cc::face_class* f = cell.faces[i];
        if(f->type != cc::INTER) continue;
        const cc::cell_class* o = cc::other_cell(f,&cell);
        if(o->index-1<0 || o->index-1>=cc::cell_num) continue; 
        const cc::vec4 dj = (o->index-1>=lo && o->index-1<hi) ? dQ[o->index-1] : dQold[o->index-1];
        double s = 2*cell.fnorm[i]-1;
        acc += cc::matpro(0.5*s*f->toface_jacobi(o->conj) + (-0.5)*f->imats.Ab + (-1)*f->imats.Tf + (-0.5*s)*f->imats.Bf, dj);
    }
    dQ[c] = cell.diag.solve4(-1*acc);
}

template <class KERNEL> inline void impl::block_sgs(KERNEL&& solve_kernel){
#pragma omp parallel
    {
#ifdef _OPENMP
        const int nt = omp_get_num_threads(), tid = omp_get_thread_num();
#else
        const int nt = 1, tid = 0;
#endif
        const int lo = (int) ((long long) cc::cell_num* tid / nt);
        const int hi = (int) ((long long) cc::cell_num* (tid+1) / nt);
        for(int c = lo; c < hi; c++) solve_kernel(c, lo, hi);
        for(int c = hi - 1; c >= lo; c--) solve_kernel(c, lo, hi); 
    }

}

static bool ifnotneedralax(cc::vec4 newed,cc::vec4 old){
    double oldein = old.e - 0.5*(old.x*old.x+old.y*old.y)/old.c;
    return (newed.c>0.5*old.c && newed.e-0.5*(newed.x*newed.x+newed.y*newed.y)/newed.c>0.5*oldein);
}

inline int impl::relax(){
    int clip = 0;
#pragma omp parallel for schedule(static) reduction(+:clip)
    allcell{
        double relaxer = 1.0;
        for(int it=0;it<6;it++){
            if(ifnotneedralax(icell(i).conser+relaxer*dQ[i], icell(i).conser)) break;
            if(it>0 && clip != 1) clip++;
            relaxer *= 0.5;
        }
        icell(i).conser = icell(i).conser + relaxer * dQ[i];
    }
    return clip;
}

inline double impl::flow_step(double dt,double cfl){
    malloc_vector();
    clear_mat();
    right_residual();
    double inv_dt = dt > 0.0 ? 1.0/dt : 0.0;
    double res = 0.0;
#pragma omp parallel for schedule(static)
    allcell form_euler_jacobi(icell(i));
#pragma omp parallel for schedule(static)
    allfac form_viscous_jacobi(iface(i));
#pragma omp parallel for schedule(static)
    allcell form_diag(icell(i), inv_dt, cfl);
#pragma omp parallel for schedule(static) reduction(+:res)
    allcell {rhs[i] = icell(i).convect - icell(i).visflux + (icell(i).vol*inv_dt)*(icell(i).conser-icell(i).lastconser);
             dQ[i].clear();res += rhs[i].c * rhs[i].c; };
    // 并行编译标记写在里面了        
    for(int sw=0;sw<urans::sweeps;sw++){Qcopy();block_sgs(solvekernel);}
    n_clip = relax();
    return sqrt(res/cc::cell_num);
}

inline void impl::turbulence_step(double dt, double cfl){
    SA::malloc_vector();
    const double inv_dt = dt > 0.0 ? 1.0 / dt : 0.0;
#pragma omp parallel for schedule(static)
    allcell SA::nurhs[i] = SA::right_residual(icell(i)) + icell(i).vol * inv_dt *
                           (icell(i).conser.c*icell(i).tur.miubl - icell(i).lastconser.c*icell(i).tur.miubl_former);
#pragma omp parallel for schedule(static)
    allcell SA::source_jacobi_diag_SA(icell(i), inv_dt, cfl);
#pragma omp parallel for schedule(static)
    allcell SA::dnu[i] = 0.0;
    for(int sw=0;sw<urans::sweeps;sw++){SA::nucopy();block_sgs(SA::solvekernel);}
    SA::relaxer();
}
