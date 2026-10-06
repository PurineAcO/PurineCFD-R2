#include "SA.hpp"
#include "classconfig.hpp"
#include "config.hpp"
#include "physic.hpp"
#include <cmath>
#include <utility>

#include "convect.hpp"
#include "dualtime.hpp"

struct mat4{
    double a[4][4] = {};
};

inline mat4 slice4(const cc::mat5& m){
    mat4 out;
    for(int r=0;r<4;r++) for(int c=0;c<4;c++) out.a[r][c] = cc::mat5_at(m,r,c);
    return out;
}

inline void add4(mat4& A,const mat4& B,const double k){
    for(int r=0;r<4;r++) for(int c=0;c<4;c++) A.a[r][c] += k*B.a[r][c];
}

inline cc::vec4 solve4(mat4 A,cc::vec4 b){
    double r[4] = {b.c,b.x,b.y,b.e};
    for(int i=0;i<4;i++){
        int p=i;
        for(int k=i+1;k<4;k++) if(std::abs(A.a[k][i])>std::abs(A.a[p][i])) p=k;
        if(p!=i){
            for(int c=0;c<4;c++) std::swap(A.a[i][c],A.a[p][c]);
            std::swap(r[i],r[p]);
        }
        if(std::abs(A.a[i][i])<1e-14) A.a[i][i] = A.a[i][i] >= 0.0 ? 1e-14 : -1e-14;
        for(int k=i+1;k<4;k++){
            const double w = A.a[k][i]/A.a[i][i];
            for(int c=i;c<4;c++) A.a[k][c] -= w*A.a[i][c];
            r[k] -= w*r[i];
        }
    }
    for(int i=3;i>=0;i--){
        for(int c=i+1;c<4;c++) r[i] -= A.a[i][c]*r[c];
        r[i] /= A.a[i][i];
    }
    return cc::vec4(r[0],r[1],r[2],r[3]);
}

inline const mat4 ident4 = []{ mat4 I; for(int r=0;r<4;r++) I.a[r][r] = 1.0; return I; }();
// T_f 为单侧两点薄层响应, 与壁面 grad_onface 的离散一致; 内部面为薄层近似
inline constexpr double k_tf = 1.0;

namespace dual_full {

inline std::vector<mat4> dmat;
inline std::vector<cc::vec4> dQold_f;

struct cell_jac{
    std::pair<cc::mat5,cc::mat5> j;
};
struct face_jac{
    mat4 Ab, Tf, Bf;
};
inline std::vector<cell_jac> cjac;
inline std::vector<face_jac> fjac;

inline double flow_step(double dt,double cfl);

// 与 dual::flow_step 同流程, 但用精确块对角 D_c 与 |A_n| 耗散
// 边界面保持标量谱半径(与 dual 一致), 内部面用精确 Jacobi
inline double flow_step(double dt,double cfl){
    const int N = cc::cell_num, F = cc::face_num;
    if((int)dmat.size()!=N){
        dmat.resize(N); dQold_f.resize(N);
        dual::rhs.resize(N); dual::dQ.resize(N); dual::dQold.resize(N);
        dual::lam.resize(N); dual::diag.resize(N);
        dual::nrhs.resize(N); dual::dnu.resize(N); dual::dnuold.resize(N);
    }
    if((int)cjac.size()!=N) cjac.resize(N);
    if((int)fjac.size()!=F) fjac.resize(F);
    dual::eval_residual();
    const double inv_dt = dt>0.0 ? 1.0/dt : 0.0;
// 单元 ∂f/∂W(按本单元面法向投影得到 A_n)、面的 |A_n|、T_f、B_f
#pragma omp parallel for schedule(static)
    for(int c=0;c<N;c++) cjac[c].j = form_pfgpw_convect(icell(c));
#pragma omp parallel for schedule(static)
    for(int i=0;i<F;i++){
        cc::face_class& f = iface(i);
        face_jac& J = fjac[i];
        const auto dm = form_pfpgq_diffusion(f);
        J.Bf = slice4(f.toface_jacobi(std::pair<cc::mat5,cc::mat5>(
            form_pfpw_from_pfpgq(dm.first,f),form_pfpw_from_pfpgq(dm.second,f))));
        if(f.type!=cc::INTER) continue;
        J.Ab = slice4(roe_abs_n(f));
        J.Tf = slice4(diffusion_grad_jac(f,dual::normal_dist(&f,*f.nei[0],f.nei[1])));
    }
    double res2 = 0.0;
#pragma omp parallel for schedule(static) reduction(+:res2)
    for(int c=0;c<N;c++){
        cc::cell_class& cell = icell(c);
        mat4& Dc = dmat[c];
        for(int r=0;r<4;r++) for(int q=0;q<4;q++) Dc.a[r][q] = 0.0;
        double lam = 0.0;
        allface(cell) lam += dual::face_lambda(cell.faces[i],cell);
        dual::lam[c] = lam;
        for(int r=0;r<4;r++) Dc.a[r][r] = cell.vol*inv_dt + lam/cfl;
        allface(cell){
            cc::face_class* f = cell.faces[i];
            const double s = 2*cell.fnorm[i]-1;
            const face_jac& J = fjac[f-&iface(0)];
            if(f->type!=cc::INTER){
                // 无粘: 标量谱半径(与 dual 一致); 粘性: 冻结梯度 Jacobi(单侧, 不乘 1/2)
                const double lamc = (std::abs(f->phy.u*f->nor.x+f->phy.v*f->nor.y)/f->len
                                      + f->otphy.a)*f->len;
                add4(Dc,ident4,0.5*lamc);
                add4(Dc,J.Bf,-s);
                // 壁面法向导数是一阶单侧差分(y1=sad), 梯度响应必须进对角
                if(f->type==cc::WALL) add4(Dc,slice4(diffusion_grad_jac(*f,cell.tur.sad)),1.0);
                continue;
            }
            add4(Dc,slice4(f->toface_jacobi(cjac[c].j)),0.5*s);
            add4(Dc,J.Ab,0.5);
            add4(Dc,J.Tf,k_tf);
            add4(Dc,J.Bf,-0.5*s);
        }
        dual::rhs[c] = cell.convect-cell.visflux;
        if(dt>0.0) dual::rhs[c] += (cell.vol*inv_dt)*(cell.conser-cell.lastconser);
        dual::dQ[c] = cc::vec4(0,0,0,0);
        res2 += dual::rhs[c].c*dual::rhs[c].c;
    }
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
            const mat4 Aj = slice4(f->toface_jacobi(cjac[j].j));
            const double arr[4] = {dj.c,dj.x,dj.y,dj.e};
            double out[4] = {0,0,0,0};
            for(int r=0;r<4;r++) for(int q=0;q<4;q++)
                out[r] += (0.5*s*Aj.a[r][q]-0.5*J.Ab.a[r][q]-k_tf*J.Tf.a[r][q]-0.5*s*J.Bf.a[r][q])*arr[q];
            acc.c += out[0];
            acc.x += out[1];
            acc.y += out[2];
            acc.e += out[3];
        }
        dual::dQ[c] = solve4(dmat[c],cc::vec4(-acc.c,-acc.x,-acc.y,-acc.e));
    };
    for(int sw=0;sw<urans::sweeps;sw++){
#pragma omp parallel for schedule(static)
        for(int c=0;c<N;c++) dQold_f[c] = dual::dQ[c];
        dual::block_sgs(relax);
    }
#pragma omp parallel for schedule(static)
    for(int c=0;c<N;c++){
        const cc::vec4 Q = icell(c).conser;
        const double ein = Q.e-0.5*(Q.x*Q.x+Q.y*Q.y)/Q.c;
        double w = 1.0;
        for(int it=0;it<6;it++){
            const cc::vec4 Qn = Q+w*dual::dQ[c];
            if(Qn.c>0.5*Q.c && Qn.e-0.5*(Qn.x*Qn.x+Qn.y*Qn.y)/Qn.c>0.5*ein) break;
            w *= 0.5;
        }
        icell(c).conser = Q+w*dual::dQ[c];
    }
    return std::sqrt(res2/N);
}

}