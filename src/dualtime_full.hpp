#include "SA.hpp"
#include "classconfig.hpp"
#include "config.hpp"
#include <cmath>
#include <cstdlib>
#include <utility>

#include "convect.hpp"
#include "dualtime.hpp"




inline constexpr double k_tf = 1.0;

namespace dual_full {

inline std::vector<cc::mat4> dmat;
inline std::vector<cc::vec4> dQold_f;
struct cell_jac{std::pair<cc::mat5,cc::mat5> j;};
struct face_jac{cc::mat4 Ab, Tf, Bf;};
inline std::vector<cell_jac> cjac;
inline std::vector<face_jac> fjac;
inline std::vector<cc::mat4> farjac;

// inline double k_dq_rel = 1.0;

// inline int guard = []{
//     const char* e = std::getenv("PURINE_GUARD");
//     return e != nullptr ? std::atoi(e) : 0;
// }();

// 上一步被欠松弛的单元数(诊断用: 看安全阀是否在持续工作)
inline int n_clip = 0;

// 一次隐式迭代
inline double flow_step(double dt,double cfl);

inline double flow_step(double dt,double cfl){
    const int N = cc::cell_num, F = cc::face_num;
    if((int)dmat.size()!=cc::cell_num){
        dmat.resize(cc::cell_num); 
        dQold_f.resize(cc::cell_num);
        dual::rhs.resize(cc::cell_num); 
        dual::dQ.resize(cc::cell_num); 
        dual::dQold.resize(cc::cell_num);
        dual::lam.resize(cc::cell_num); 
        dual::diag.resize(cc::cell_num);
        dual::nrhs.resize(cc::cell_num); 
        dual::dnu.resize(cc::cell_num); 
        dual::dnuold.resize(cc::cell_num);
    }
    if((int)cjac.size()!=cc::cell_num) cjac.resize(cc::cell_num);
    if((int)fjac.size()!=cc::face_num) fjac.resize(cc::face_num);
    if((int)farjac.size()!=cc::cell_num) farjac.resize(cc::cell_num);

    // 形成右侧残差
    dual::eval_residual();
    const double inv_dt = dt>0.0 ? 1.0/dt : 0.0;

#pragma omp parallel for schedule(static)
    // 首先形成各个格子的无粘通量Jacobi
    for(int c=0;c<cc::cell_num;c++){
        cjac[c].j = form_pfgpw_convect(icell(c));
        for(int r=0;r<4;r++) for(int q=0;q<4;q++) farjac[c].a[r][q] = 0.0;
        cc::cell_class& cell = icell(c);
        allface(cell){
            cc::face_class* f = cell.faces[i];
            if(f->type != cc::FAR) continue;
            farjac[c] += far_convect_jacobi(f);
        }
    }

#pragma omp parallel for schedule(static)
    // 形成各个面上的Roe Jacobi和两组SA Jacobi
    for(int i=0;i<F;i++){
        cc::face_class& f = iface(i);
        face_jac& J = fjac[i];
        const auto dm = SA::form_pfpgq_diffusion(f);
        J.Bf = cc::slice4(f.toface_jacobi(std::pair<cc::mat5,cc::mat5>(
            SA::form_pfpw_from_pfpgq(dm.first,f),SA::form_pfpw_from_pfpgq(dm.second,f))));
        if(f.type!=cc::INTER) continue;
        J.Ab = cc::slice4(roe_abs_n(f));
        J.Tf = cc::slice4(SA::diffusion_grad_jac(f,dual::normal_dist(&f,*f.nei[0],f.nei[1])));
    }
    double res2 = 0.0;  // 残差,通过连续性方程度量
    
#pragma omp parallel for schedule(static) reduction(+:res2)
    // 组装Jacobi大矩阵对角区域
    allcell{
        cc::cell_class& cell = icell(i);
        cc::mat4& Dc = dmat[i];
        for(int r=0;r<4;r++) for(int q=0;q<4;q++) Dc.clear();
        double lam = 0.0;
        allface(cell) lam += dual::face_lambda(cell.faces[i],cell);
        dual::lam[i] = lam;
        Dc += (cell.vol*inv_dt + lam/cfl) * cc::eyes;
        allface(cell){
            cc::face_class* f = cell.faces[i];
            const double s = 2*cell.fnorm[i]-1;
            const face_jac& J = fjac[f-&iface(0)];
            if(f->type!=cc::INTER){
                if(f->type==cc::WALL){
                    const double lamc = (std::abs(f->phy.u*f->nor.x+f->phy.v*f->nor.y)/f->len
                                          + f->otphy.a)*f->len;
                    Dc += 0.5*lamc*cc::eyes;
                    Dc += cc::slice4(SA::diffusion_grad_jac(*f,cell.tur.sad));
                }
                Dc += (-s)*J.Bf;
                continue;
            }
            Dc += cc::slice4(0.5*s*f->toface_jacobi(cjac[i].j)) + 0.5*J.Ab + k_tf*J.Tf + (-0.5*s)*J.Bf;
        }
        Dc += farjac[i];
        dual::rhs[i] = cell.convect-cell.visflux;
        if(dt>0.0) dual::rhs[i] += (cell.vol*inv_dt)*(cell.conser-cell.lastconser);
        dual::dQ[i] = cc::vec4(0,0,0,0);
        res2 += dual::rhs[i].c*dual::rhs[i].c;
    }

    // 按块分配求解核心,进行一次GS
    auto relax = [](int c,int lo,int hi){
        const cc::cell_class& cell = icell(c);
        cc::vec4 acc = dual::rhs[c];
        allface(cell){
            const cc::face_class* f = cell.faces[i];
            if(f->type!=cc::INTER) continue;
            const int j = dual::other_cell(f,cell)->index-1;
            if(j<0 || j>=cc::cell_num) continue;
            const cc::vec4 dj = (j>=lo && j<hi) ? dual::dQ[j] : dQold_f[j];
            const double s = 2*cell.fnorm[i]-1;
            const face_jac& J = fjac[f-&iface(0)];
            const cc::vec4 arr = {dj.c,dj.x,dj.y,dj.e};
            cc::vec4 out = cc::matpro(0.5*s*slice4(f->toface_jacobi(cjac[j].j))+(-0.5)*J.Ab+(-k_tf)*J.Tf+(-0.5*s)*J.Bf,arr);
            acc += out;
        }
        dual::dQ[c] = dmat[c].solve4(-1*acc);
    };

    // 重复扫描
    for(int sw=0;sw<urans::sweeps;sw++){
#pragma omp parallel for schedule(static)
        for(int c=0;c<N;c++) dQold_f[c] = dual::dQ[c];
        dual::block_sgs(relax);
    }

    // 更新守恒量,同时进行松弛
    int clip = 0;
#pragma omp parallel for schedule(static) reduction(+:clip)
    for(int c=0;c<N;c++){
        const cc::vec4 Q = icell(c).conser;
        const double ein = Q.e-0.5*(Q.x*Q.x+Q.y*Q.y)/Q.c;
        cc::vec4 dQ = dual::dQ[c];
        // if((guard&2) && (!std::isfinite(dQ.c) || !std::isfinite(dQ.x) ||
        //                  !std::isfinite(dQ.y) || !std::isfinite(dQ.e))){
        //     dQ = cc::vec4(0,0,0,0);
        //     ++clip;
        // }
        double w = 1.0;
        // if(guard&1){
        //     const double qv[4] = {Q.c,Q.x,Q.y,Q.e};
        //     const double dv[4] = {dQ.c,dQ.x,dQ.y,dQ.e};
        //     double mx = 0.0;
        //     for(int k=0;k<4;k++) mx = std::max(mx,std::abs(dv[k])/std::max(std::abs(qv[k]),1e-30));
        //     if(mx > k_dq_rel){ w = k_dq_rel/mx; ++clip; }
        // }
        for(int it=0;it<6;it++){
            const cc::vec4 Qn = Q+w*dQ;
            if(Qn.c>0.5*Q.c && Qn.e-0.5*(Qn.x*Qn.x+Qn.y*Qn.y)/Qn.c>0.5*ein) break;
            w *= 0.5;
        }
        icell(c).conser = Q+w*dQ;
    }
    n_clip = clip;
    return std::sqrt(res2/N);
}

}