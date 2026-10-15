#include "SA.hpp"
#include "classconfig.hpp"
#include "config.hpp"
#include <cmath>
#include <cstdlib>
#include <utility>

#include "convect.hpp"
#include "dualtime.hpp"


inline cc::mat4 ident4 = []{ cc::mat4 I; for(int r=0;r<4;r++) I.a[r][r] = 1.0; return I; }();

inline constexpr double k_tf = 1.0;

namespace dual_full {

inline std::vector<cc::mat4> dmat;
inline std::vector<cc::vec4> dQold_f;
struct cell_jac{std::pair<cc::mat5,cc::mat5> j;};
struct face_jac{cc::mat4 Ab, Tf, Bf;};
inline std::vector<cell_jac> cjac;
inline std::vector<face_jac> fjac;
inline std::vector<cc::mat4> farjac;

// 远场面对角块: true 用数值微分的 ∂R/∂Q_c, false 退回标量谱半径 0.5*lam*I
// 环境变量 PURINE_FARJAC=0 可关闭(对比实验用)
// inline bool far_jac = []{
//     const char* e = std::getenv("PURINE_FARJAC");
//     return !(e != nullptr && e[0] == '0');
// }();

// 单次内迭代各守恒分量的最大相对变化: 残差尖峰时 dQ 会跳出量级,
// 超出就按比例欠松弛(等效自适应 ω), 防止单点内迭代被一步打崩
// 阀值要远高于正常工作点(定常/内迭代一般是 10%~50%), 所以取 100% 只拦异常
inline double k_dq_rel = 1.0;

// 更新保护开关位掩码, 默认 0(全关 = 原始行为), 用 PURINE_GUARD 打开:
//   bit0 守恒量更新的幅值限幅(自适应欠松弛)
//   bit1 dQ 非有限时丢弃该单元的更新
inline int guard = []{
    const char* e = std::getenv("PURINE_GUARD");
    return e != nullptr ? std::atoi(e) : 0;
}();

// 上一步被欠松弛的单元数(诊断用: 看安全阀是否在持续工作)
inline int n_clip = 0;

// 远场面的无粘通量是 euler_flux(Q_bc(Q_c), n), Q_bc 由内点的特征关系给出,
// 故 ∂R_c/∂Q_c 是满 4x4, 标量谱半径只是它的替身。这里用中心差分直接量出来。
// inline double q_at(const cc::vec4& Q,int k){
//     switch(k){
//         case 0: return Q.c;
//         case 1: return Q.x;
//         case 2: return Q.y;
//         default: return Q.e;
//     }
// }

// inline void q_set(cc::vec4& Q,int k,double v){
//     switch(k){
//         case 0: Q.c = v; break;
//         case 1: Q.x = v; break;
//         case 2: Q.y = v; break;
//         default: Q.e = v; break;
//     }
// }

// // 远场面的 BC 面值, 与 boundary.hpp 的 far_field_boundary 逐行一致
// inline cc::vecp far_bc_state(const cc::cell_class& c,const cc::face_class& f){
//     const double rho_inf = FAR_DEFINE.p/(cc::R*FAR_DEFINE.T);
//     const double a_inf = get_sonic_velocity(FAR_DEFINE.T);
//     double nx = f.nor.x/f.len, ny = f.nor.y/f.len;
//     if(nx*(f.mid.x-c.center.x)+ny*(f.mid.y-c.center.y)<0.0){ nx=-nx; ny=-ny; }
//     const double a = get_sonic_velocity(c.phy.T);
//     const double vn = c.phy.u*nx+c.phy.v*ny;
//     const double vt = -c.phy.u*ny+c.phy.v*nx;
//     const double vn_inf = FAR_DEFINE.u*nx+FAR_DEFINE.v*ny;
//     const double vt_inf = -FAR_DEFINE.u*ny+FAR_DEFINE.v*nx;
//     const double Rp = vn+2.0*a/(cc::gamma-1.0);
//     const double Rm = vn_inf-2.0*a_inf/(cc::gamma-1.0);
//     const double vn_star = 0.5*(Rp+Rm);
//     const double a_star = 0.25*(cc::gamma-1.0)*(Rp-Rm);
//     double s, vt_star;
//     if(vn_star>=0.0){ s = c.otphy.p/std::pow(c.phy.rho,cc::gamma); vt_star = vt; }
//     else{ s = FAR_DEFINE.p/std::pow(rho_inf,cc::gamma); vt_star = vt_inf; }
//     const double rho_b = std::pow(a_star*a_star/(cc::gamma*s),1.0/(cc::gamma-1.0));
//     const double p_b = s*std::pow(rho_b,cc::gamma);
//     return cc::vecp(rho_b,vn_star*nx-vt_star*ny,vn_star*ny+vt_star*nx,p_b/(cc::R*rho_b));
// }

// // 远场面的无粘通量(带面长), 与 convect_ROE 在 L=R 时的退化结果逐位一致
// // 注意必须沿用 otphy 的 p 与 e, 不能从守恒量反算: Cp-Cv 与 R 并不严格相等
// inline cc::vec4 far_face_flux(const cc::face_class& f,const cc::vecp& bc){
//     const double sgn = f.nei[0]!=nullptr
//         ? (cc::dot(f.nor,f.mid-f.nei[0]->center)>0 ? 1.0 : -1.0)
//         : (cc::dot(f.nor,f.nei[1]->center-f.mid)>0 ? 1.0 : -1.0);
//     const double nx = f.nor.x*sgn, ny = f.nor.y*sgn;
//     const double p = bc.rho*cc::R*bc.T;
//     const double e = cc::Cv*bc.T+0.5*(bc.u*bc.u+bc.v*bc.v);
//     const double un = bc.u*nx+bc.v*ny;
//     return cc::vec4(bc.rho*un,bc.rho*bc.u*un+p*nx,bc.rho*bc.v*un+p*ny,un*(bc.rho*e+p));
// }

// // 单元所有远场面的净流出无粘通量之和, 符号与 assemble_flux 一致
// inline cc::vec4 far_flux_sum(cc::cell_class& cell){
//     cc::vec4 sum(0,0,0,0);
//     allface(cell){
//         cc::face_class* f = cell.faces[i];
//         if(f->type==cc::INTER || f->type==cc::WALL) continue;
//         sum += (f->nei[0]==&cell ? 1.0 : -1.0)*far_face_flux(*f,far_bc_state(cell,*f));
//     }
//     return sum;
// }

// // 远场面对角块 Jacobian = 中心差分 ∂(Σ far_face_flux)/∂Q_c; 该单元无远场面时 J 保持 0
// inline void far_face_jacobian(cc::cell_class& cell,mat4& J,double rel = 1e-6){
//     bool has_far = false;
//     allface(cell){
//         const cc::face_class* f = cell.faces[i];
//         if(f->type!=cc::INTER && f->type!=cc::WALL){ has_far = true; break; }
//     }
//     if(!has_far) return;
//     const cc::vec4 base = cell.conser;
//     for(int q=0;q<4;q++){
//         double eps = rel*std::max(std::abs(q_at(base,q)),1e-30);
//         bool ok = false;
//         // 扰动过大可能让扰动后的状态失效(T<=0 -> a 为 nan), 减半步长重试;
//         // 仍失败则该列置 0, 宁缺勿把 nan 带进对角块
//         for(int t=0;t<4 && !ok;t++){
//             cc::vec4 Qp = base, Qm = base;
//             q_set(Qp,q,q_at(base,q)+eps);
//             q_set(Qm,q,q_at(base,q)-eps);
//             cell.conser = Qp; cell.prim();
//             const bool okp = std::isfinite(cell.otphy.a) && cell.phy.T > 0.0;
//             const cc::vec4 Rp = okp ? far_flux_sum(cell) : cc::vec4(0,0,0,0);
//             cell.conser = Qm; cell.prim();
//             const bool okm = std::isfinite(cell.otphy.a) && cell.phy.T > 0.0;
//             const cc::vec4 Rm = okm ? far_flux_sum(cell) : cc::vec4(0,0,0,0);
//             if(okp && okm){
//                 const double rp[4] = {Rp.c,Rp.x,Rp.y,Rp.e};
//                 const double rm[4] = {Rm.c,Rm.x,Rm.y,Rm.e};
//                 for(int r=0;r<4;r++) J.a[r][q] = (rp[r]-rm[r])/(2.0*eps);
//                 ok = true;
//             }else{
//                 eps *= 0.5;
//             }
//         }
//         if(!ok) for(int r=0;r<4;r++) J.a[r][q] = 0.0;
//     }
//     cell.conser = base;
//     cell.prim();
// }

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
        // for(int r=0;r<4;r++) for(int q=0;q<4;q++) farjac[c].a[r][q] = 0.0;
        // if(far_jac) far_face_jacobian(icell(c),farjac[c]);
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
        for(int r=0;r<4;r++) Dc.a[r][r] = cell.vol*inv_dt + lam/cfl;
        allface(cell){
            cc::face_class* f = cell.faces[i];
            const double s = 2*cell.fnorm[i]-1;
            const face_jac& J = fjac[f-&iface(0)];
            // if(f->type!=cc::INTER){
            //     // 壁面无粘: 标量谱半径; 远场: 由 far_face_jacobian 统一加(见循环后)
            //     if(f->type==cc::WALL){
            //         const double lamc = (std::abs(f->phy.u*f->nor.x+f->phy.v*f->nor.y)/f->len
            //                               + f->otphy.a)*f->len;
            //         add4(Dc,ident4,0.5*lamc);
            //         // 壁面法向导数是一阶单侧差分(y1=sad), 梯度响应必须进对角
            //         add4(Dc,slice4(SA::diffusion_grad_jac(*f,cell.tur.sad)),1.0);
            //     }else if(!far_jac){
            //         const double lamc = (std::abs(f->phy.u*f->nor.x+f->phy.v*f->nor.y)/f->len
            //                               + f->otphy.a)*f->len;
            //         add4(Dc,ident4,0.5*lamc);
            //     }
            //     // 粘性: 冻结梯度 Jacobi(单侧, 不乘 1/2)
            //     add4(Dc,J.Bf,-s);
            //     continue;
            // }

            Dc += cc::slice4(0.5*s*f->toface_jacobi(cjac[i].j)) + 0.5*J.Ab + k_tf*J.Tf + (-0.5*s)*J.Bf;
            // add4(Dc,slice4(f->toface_jacobi(cjac[i].j)),0.5*s);
            // add4(Dc,J.Ab,0.5);
            // add4(Dc,J.Tf,k_tf);
            // add4(Dc,J.Bf,-0.5*s);
        }
        // 远场面的 ∂R/∂Q_c 是满 4x4, 且天然只依赖本单元, 整块加进对角
        // if(far_jac) add4(Dc,farjac[c],1.0);
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
        if((guard&2) && (!std::isfinite(dQ.c) || !std::isfinite(dQ.x) ||
                         !std::isfinite(dQ.y) || !std::isfinite(dQ.e))){
            dQ = cc::vec4(0,0,0,0);
            ++clip;
        }
        double w = 1.0;
        if(guard&1){
            const double qv[4] = {Q.c,Q.x,Q.y,Q.e};
            const double dv[4] = {dQ.c,dQ.x,dQ.y,dQ.e};
            double mx = 0.0;
            for(int k=0;k<4;k++) mx = std::max(mx,std::abs(dv[k])/std::max(std::abs(qv[k]),1e-30));
            if(mx > k_dq_rel){ w = k_dq_rel/mx; ++clip; }
        }
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