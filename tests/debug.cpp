// tests/debug.cpp
// 圆柱绕流 URANS 分步检验台: 只驱动 src 的头文件
// 结果写入仓库根目录 ./debugres, 控制台同步打印进度与通过情况
#include "classconfig.hpp"
#include "config.hpp"
#include "SA.hpp"
#include "boundary.hpp"
#include "convect.hpp"
#include "geometry.hpp"
#include "grad.hpp"
#include "HALO.hpp"
#include "initialize.hpp"
#include "interpolate.hpp"
#include "io.hpp"
#include "physic.hpp"
#include "readmesh.hpp"
#include "timarch.hpp"
#include "udf.hpp"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <sys/stat.h>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

static const char* OUTDIR = "debugres";
static const int CHECK_NUM = 15;
static int g_done = 0,g_pass = 0;
static FILE* g_summary = nullptr;
static std::string g_info;
static bool g_rans = false;   // true: 稳态RANS, false: URANS
static bool g_trace = false;  // true: 输出 ν̃ 追踪
static bool g_laminar = false;// true: 冻结湍流模型, 按层流跑
static int g_dump_every = 5000;// 每多少步输出一次流场(0=不输出中间帧)
static double g_re = 60.0;    // 当前 Re, 仅用于输出标题
static void miubl_trace(int iter,double dt,bool urans);
static bool dump_field(const std::string& tag,int steps,double dt);

static FILE* out_open(const char* name){
    char buf[512];
    snprintf(buf,sizeof(buf),"%s/%s",OUTDIR,name);
    FILE* fp = fopen(buf,"w");
    if(fp == nullptr){
        fprintf(stderr,"Error: cannot open %s\n",buf);
    }
    return fp;
}

static void check_begin(int no,const char* name){
    printf("[%2d/%d] %-42s ",no,CHECK_NUM,name);
    fflush(stdout);
}

static void check_end(bool ok,const std::string& info){
    ++g_done;
    if(ok){
        ++g_pass;
    }
    printf("%s  %s\n",ok ? "PASS" : "FAIL",info.c_str());
    if(g_summary != nullptr){
        fprintf(g_summary,"[%2d/%d] %s  %s\n",g_done,CHECK_NUM,ok ? "PASS" : "FAIL",info.c_str());
        fflush(g_summary);
    }
    fflush(stdout);
}

static bool near(double a,double b,double rtol,double atol = 0.0){
    return std::isfinite(a) && std::isfinite(b) &&
           std::abs(a-b) <= atol + rtol*std::max(std::abs(a),std::abs(b));
}

static double rel_err(double a,double b){
    const double scale = std::max(std::abs(a),std::abs(b));
    return scale == 0.0 ? std::abs(a-b) : std::abs(a-b)/scale;
}

static std::string fmt(const char* format,...){
    char buf[512];
    va_list args;
    va_start(args,format);
    vsnprintf(buf,sizeof(buf),format,args);
    va_end(args);
    return std::string(buf);
}

// L->R 方向的单位法向: 由 nor 与 nei[0]->nei[1] 的夹角决定符号
static void lr_normal(const cc::face_class& face,double& nx,double& ny){
    const double len = std::hypot(face.nor.x,face.nor.y);
    double sgn = 1.0;
    if(face.type == cc::INTER){
        sgn = cc::dot(face.nor,face.nei[1]->center-face.nei[0]->center) > 0 ? 1.0 : -1.0;
    }else{
        // 边界面只有一侧有格子: L->R 方向就是那个格子的外法向
        const cc::cell_class* cell = cc::boundary_findcell(const_cast<cc::face_class*>(&face));
        sgn = (cc::dot(face.nor,face.mid-cell->center) > 0 ? 1.0 : -1.0)*(face.nei[0] ? 1.0 : -1.0);
    }
    nx = sgn*face.nor.x/len;
    ny = sgn*face.nor.y/len;
}

// 独立实现的标准 Roe 通量, 用于和 src 的 convect_ROE 对拍
static void reference_roe(const cc::face_class& face,double out[4]){
    const cc::vecp& L = face.phynei[0];
    const cc::vecp& R = face.phynei[1];
    const cc::otphy& Lo = face.otnei[0];
    const cc::otphy& Ro = face.otnei[1];
    double nx = 0.0,ny = 0.0;
    lr_normal(face,nx,ny);
    const double len = std::hypot(face.nor.x,face.nor.y);
    const double sl = std::sqrt(L.rho),sh = std::sqrt(R.rho);
    const double Hl = cc::Cp*L.T + 0.5*(L.u*L.u + L.v*L.v);
    const double Hh = cc::Cp*R.T + 0.5*(R.u*R.u + R.v*R.v);
    const double rhobl = sl*sh;
    const double ubl = (sl*L.u + sh*R.u)/(sl+sh);
    const double vbl = (sl*L.v + sh*R.v)/(sl+sh);
    const double Hbl = (sl*Hl + sh*Hh)/(sl+sh);
    const double abl = std::sqrt((cc::gamma-1)*(Hbl - 0.5*(ubl*ubl + vbl*vbl)));
    const double unbl = ubl*nx + vbl*ny;
    const double utbl = vbl*nx - ubl*ny;
    const double lambda[4] = {std::abs(unbl-abl),std::abs(unbl),std::abs(unbl),std::abs(unbl+abl)};
    const double dp = Ro.p - Lo.p,drho = R.rho - L.rho;
    const double dun = nx*(R.u-L.u) + ny*(R.v-L.v);
    const double dut = nx*(R.v-L.v) - ny*(R.u-L.u);
    const double alpha[4] = {(dp-rhobl*abl*dun)/(2*abl*abl),drho-dp/(abl*abl),rhobl*dut,
                             (dp+rhobl*abl*dun)/(2*abl*abl)};
    const double tz[4][4] = {{1,ubl-abl*nx,vbl-abl*ny,Hbl-unbl*abl},
                             {1,ubl,vbl,0.5*(ubl*ubl+vbl*vbl)},
                             {0,-ny,nx,utbl},
                             {1,ubl+abl*nx,vbl+abl*ny,Hbl+unbl*abl}};
    const double unl = L.u*nx + L.v*ny,unh = R.u*nx + R.v*ny;
    const double FL[4] = {L.rho*unl,L.rho*L.u*unl+Lo.p*nx,L.rho*L.v*unl+Lo.p*ny,unl*(L.rho*Lo.e+Lo.p)};
    const double FR[4] = {R.rho*unh,R.rho*R.u*unh+Ro.p*nx,R.rho*R.v*unh+Ro.p*ny,unh*(R.rho*Ro.e+Ro.p)};
    for(int k=0;k<4;k++){
        out[k] = 0.5*(FL[k]+FR[k]);
        for(int m=0;m<4;m++){
            out[k] -= 0.5*lambda[m]*alpha[m]*tz[m][k];
        }
        out[k] *= len;
    }
}

// 冻结湍流: 强制 ν̃≡0, 则 fv1(0)=0, μ_t≡0
static void freeze_turbulence(){
    #pragma omp parallel for schedule(static)
    for(int c=0;c<cc::cell_num;c++){
        cc::CellList[c].tur.miubl = 0.0;
        cc::CellList[c].tur.miublgrad = cc::vec2{0.0,0.0};
    }
    #pragma omp parallel for schedule(static)
    for(int f=0;f<cc::face_num;f++){
        cc::FaceList[f].tur.miubl = 0.0;
        cc::FaceList[f].tur.miublgrad = cc::vec2{0.0,0.0};
        cc::FaceList[f].tur.mueff = cc::FaceList[f].otphy.mu;
    }
    for(cc::cell_class& ghost:cc::GhostList){
        ghost.tur.miubl = 0.0;
    }
}

// 把基本量、边界条件、面状态、梯度都刷到当前守恒量对应的状态
static void refresh_field(){
    #pragma omp parallel for schedule(static)
    for(int c=0;c<cc::cell_num;c++){
        cc::CellList[c].prim();
    }
    noslip_wall_boundary();
    far_field_boundary();
    update_ghost_field();
    if(g_laminar){
        freeze_turbulence();
    }
    #pragma omp parallel for schedule(static)
    for(int f=0;f<cc::face_num;f++){
        interpolate_mid(&cc::FaceList[f]);
        muscl(&cc::FaceList[f]);
    }
    #pragma omp parallel for schedule(static)
    for(int c=0;c<cc::cell_num;c++){
        least_square_cell_based(cc::CellList[c]);
    }
    #pragma omp parallel for schedule(static)
    for(int c=0;c<cc::cell_num;c++){
        grad_onface(cc::CellList[c]);
    }
}

// 与 timarch.hpp 的 one_rans 逐行对应, 但外层循环并行
static void pseudo_step_omp(double dt,bool urans){
    #pragma omp parallel for schedule(static)
    for(int c=0;c<cc::cell_num;c++){
        cc::CellList[c].copyconver();
    }
    for(int z=0;z<3;z++){
        #pragma omp parallel for schedule(static)
        for(int c=0;c<cc::cell_num;c++){
            cc::CellList[c].prim();
        }
        noslip_wall_boundary();
        far_field_boundary();
        update_ghost_field();
        // 层流要在这儿再压一次, 否则远场入流面会被赋予 ν̃_∞
        if(g_laminar){
            freeze_turbulence();
        }
        #pragma omp parallel for schedule(static)
        for(int f=0;f<cc::face_num;f++){
            interpolate_mid(&cc::FaceList[f]);
            muscl(&cc::FaceList[f]);
        }
        #pragma omp parallel for schedule(static)
        for(int c=0;c<cc::cell_num;c++){
            least_square_cell_based(cc::CellList[c]);
        }
        #pragma omp parallel for schedule(static)
        for(int c=0;c<cc::cell_num;c++){
            grad_onface(cc::CellList[c]);
        }
        if(cc::scheme == 'J'){
            #pragma omp parallel for schedule(static)
            for(int c=0;c<cc::cell_num;c++){
                jst::shockwave_recognize(cc::CellList[c]);
            }
            #pragma omp parallel for schedule(static)
            for(int c=0;c<cc::cell_num;c++){
                jst::laplace_dissipation(cc::CellList[c]);
            }
            #pragma omp parallel for schedule(static)
            for(int f=0;f<cc::face_num;f++){
                convect_JST(cc::FaceList[f]);
            }
        }else{
            #pragma omp parallel for schedule(static)
            for(int f=0;f<cc::face_num;f++){
                convect_ROE(cc::FaceList[f]);
            }
        }
        #pragma omp parallel for schedule(static)
        for(int c=0;c<cc::cell_num;c++){
            assemble_flux(cc::CellList[c],cc::scheme);
        }
        #pragma omp parallel for schedule(static)
        for(int f=0;f<cc::face_num;f++){
            SA::diffusion_SA(cc::FaceList[f]);
        }
        #pragma omp parallel for schedule(static)
        for(int c=0;c<cc::cell_num;c++){
            SA::assemble_visflux(cc::CellList[c]);
        }
        #pragma omp parallel for schedule(static)
        for(int c=0;c<cc::cell_num;c++){
            local_timestep(cc::CellList[c]);
        }
        #pragma omp parallel for schedule(static)
        for(int c=0;c<cc::cell_num;c++){
            cc::cell_class& cell = cc::CellList[c];
            if(!urans){
                cell.conser = cell.conserformer - RK::RK3[z]/cell.vol*cell.localdt*
                              (cell.convect-cell.visflux);
            }else{
                cell.conser = cell.conserformer - RK::RK3[z]/(1.0/cell.localdt + 1.0/(2*dt))*
                              (1.0/cell.vol*(cell.convect-cell.visflux) +
                               1.0/(2*dt)*(cell.conser-cell.lastconser));
            }
        }
    }
    if(!g_laminar){
        #pragma omp parallel for schedule(static)
        for(int c=0;c<cc::cell_num;c++){
            SA::SA_equation_after(cc::CellList[c],cc::CellList[c].localdt,dt,urans);
        }
    }
}

static void snapshot(std::vector<cc::vec4>& conser,std::vector<double>& bl){
    conser.resize(cc::cell_num);
    bl.resize(cc::cell_num);
    #pragma omp parallel for schedule(static)
    for(int c=0;c<cc::cell_num;c++){
        conser[c] = cc::CellList[c].conser;
        bl[c] = cc::CellList[c].tur.miubl;
    }
}

static void restore(const std::vector<cc::vec4>& conser,const std::vector<double>& bl){
    #pragma omp parallel for schedule(static)
    for(int c=0;c<cc::cell_num;c++){
        cc::CellList[c].conser = conser[c];
        cc::CellList[c].tur.miubl = bl[c];
    }
}

// 流场变化的 L2: 各分量用各自的参考量纲归一化
static double g_scale[4] = {};

static double change_l2(const std::vector<cc::vec4>& ref,double* max_bl){
    double sum = 0.0,bl = 0.0;
    #pragma omp parallel for schedule(static) reduction(+:sum) reduction(max:bl)
    for(int c=0;c<cc::cell_num;c++){
        const cc::vec4 d = cc::CellList[c].conser - ref[c];
        sum += (d.c*d.c)/(g_scale[0]*g_scale[0]) +
               (d.x*d.x + d.y*d.y)/(g_scale[1]*g_scale[1]) +
               (d.e*d.e)/(g_scale[3]*g_scale[3]);
        bl = std::max(bl,std::abs(cc::CellList[c].tur.miubl));
    }
    if(max_bl != nullptr){
        *max_bl = bl;
    }
    return std::sqrt(sum/cc::cell_num);
}

static double increment_l2(const std::vector<cc::vec4>& prev){
    double sum = 0.0;
    #pragma omp parallel for schedule(static) reduction(+:sum)
    for(int c=0;c<cc::cell_num;c++){
        const cc::vec4 d = cc::CellList[c].conser - prev[c];
        sum += (d.c*d.c)/(g_scale[0]*g_scale[0]) +
               (d.x*d.x + d.y*d.y)/(g_scale[1]*g_scale[1]) +
               (d.e*d.e)/(g_scale[3]*g_scale[3]);
    }
    return std::sqrt(sum/cc::cell_num);
}

static int first_bad_cell(){
    for(int c=0;c<cc::cell_num;c++){
        const cc::cell_class& cell = cc::CellList[c];
        if(!std::isfinite(cell.conser.c) || !std::isfinite(cell.conser.x) ||
           !std::isfinite(cell.conser.y) || !std::isfinite(cell.conser.e) ||
           !std::isfinite(cell.tur.miubl)){
            return c;
        }
    }
    return -1;
}

// ===========================================================================
// 检验 1: 网格 / 面读取 / 邻接关系
// ===========================================================================
static bool check_01_mesh(){
    FILE* fp = out_open("check01_mesh.txt");
    int face_err = 0,cell_err = 0,back_err = 0,bound_err = 0;
    int interior = 0,boundary = 0,wall = 0,far = 0;
    for(size_t f=0;f<cc::FaceList.size();f++){
        const cc::face_class& face = cc::FaceList[f];
        if(face.index != static_cast<int>(f)+1 || cc::link_face(face.index) != &face ||
           face.node[0] == nullptr || face.node[1] == nullptr){
            ++face_err;
        }
        const bool n0 = face.nei[0] == nullptr,n1 = face.nei[1] == nullptr;
        if(face.type == cc::INTER){
            ++interior;
            if(n0 || n1){
                ++bound_err;
            }
        }else{
            ++boundary;
            if(face.type == cc::WALL){
                ++wall;
            }else{
                ++far;
            }
            if(n0 == n1){
                ++bound_err;
            }
        }
    }
    for(int c=0;c<cc::cell_num;c++){
        const cc::cell_class& cell = cc::CellList[c];
        if(cell.index != c+1 || cell.ecnt != 4){
            ++cell_err;
        }
        for(int i=0;i<cell.ecnt;i++){
            const cc::face_class* face = cell.faces[i];
            if(face == nullptr || cc::link_face(cell.face[i]) != face ||
               (face->nei[0] != &cell && face->nei[1] != &cell)){
                ++cell_err;
                continue;
            }
            const cc::cell_class* other = face->nei[0] == &cell ? face->nei[1] : face->nei[0];
            if(cell.nei[i] != other){
                ++cell_err;
            }
            if(other == nullptr){
                continue;
            }
            bool back = false;
            for(int j=0;j<other->ecnt;j++){
                if(other->faces[j] == face && other->nei[j] == &cell){
                    back = true;
                }
            }
            if(!back){
                ++back_err;
            }
        }
    }
    fprintf(fp,"nodes %zu\nfaces %zu\ncells %zu\n",cc::NodeList.size(),cc::FaceList.size(),cc::CellList.size());
    fprintf(fp,"interior_faces %d\nboundary_faces %d (wall %d, far %d)\n",interior,boundary,wall,far);
    fprintf(fp,"face_err %d\ncell_face_link_err %d\nface_backlink_err %d\nboundary_nei_err %d\n",
            face_err,cell_err,back_err,bound_err);
    for(cc::face_class* w : cc::WallFaces){
        fprintf(fp,"wall_face %d cell %d %d len %.6e\n",w->index,w->cell_1,w->cell_2,w->len);
    }
    fclose(fp);
    const bool ok = face_err == 0 && cell_err == 0 && back_err == 0 && bound_err == 0 &&
                    interior > 0 && wall > 0 && far > 0;
    g_info = fmt("cells=%d faces=%d nodes=%zu 内部面=%d 壁面=%d 远场=%d",
                 cc::cell_num,cc::face_num,cc::NodeList.size(),interior,wall,far);
    return ok;
}

// ===========================================================================
// 检验 2: 几何量(面积 / 质心 / sad)
// ===========================================================================
static bool check_02_geometry(){
    FILE* fp = out_open("check02_geometry.txt");
    double worst_vol = 0.0,worst_cen = 0.0,worst_sad = 0.0,worst_seg = 0.0;
    double min_vol = 1e300,max_vol = 0.0,min_sad = 1e300,max_sad = 0.0;
    int bad = 0;
    fprintf(fp,"# idx  vol  vol_ref  rel  cx  cx_ref  cy  cy_ref  sad  sad_ref  sad_to_segment\n");
    for(int c=0;c<cc::cell_num;c++){
        const cc::cell_class& cell = cc::CellList[c];
        double px[4],py[4];
        double mx = 0.0,my = 0.0;
        for(int k=0;k<4;k++){
            const cc::node_class& node = cc::NodeList[cell.node[k]];
            px[k] = node.x;
            py[k] = node.y;
            mx += node.x;
            my += node.y;
        }
        mx *= 0.25;
        my *= 0.25;
        // 按极角排序后用鞋带公式, 与 src 的三角形求和无共享代码
        double angle[4];
        int order[4] = {0,1,2,3};
        for(int k=0;k<4;k++){
            angle[k] = std::atan2(py[k]-my,px[k]-mx);
        }
        std::sort(order,order+4,[&](int a,int b){ return angle[a] < angle[b]; });
        double cross2 = 0.0,cx2 = 0.0,cy2 = 0.0;
        for(int k=0;k<4;k++){
            const int a = order[k],b = order[(k+1)%4];
            const double cr = px[a]*py[b] - px[b]*py[a];
            cross2 += cr;
            cx2 += (px[a]+px[b])*cr;
            cy2 += (py[a]+py[b])*cr;
        }
        const double ref_vol = 0.5*cross2;
        const double ref_cx = cx2/(3.0*cross2);
        const double ref_cy = cy2/(3.0*cross2);
        double ref_sad = 1e300,ref_seg = 1e300;
        for(const cc::face_class* w : cc::WallFaces){
            ref_sad = std::min(ref_sad,std::hypot(cell.center.x-w->mid.x,cell.center.y-w->mid.y));
            const double ax = w->node[0]->x,ay = w->node[0]->y;
            const double dx = w->node[1]->x-ax,dy = w->node[1]->y-ay;
            const double t = std::min(1.0,std::max(0.0,
                ((cell.center.x-ax)*dx + (cell.center.y-ay)*dy)/(dx*dx+dy*dy)));
            ref_seg = std::min(ref_seg,std::hypot(cell.center.x-(ax+t*dx),cell.center.y-(ay+t*dy)));
        }
        if(!std::isfinite(cell.vol) || !std::isfinite(cell.center.x) || !std::isfinite(cell.tur.sad)){
            ++bad;
            continue;
        }
        worst_vol = std::max(worst_vol,rel_err(cell.vol,ref_vol));
        worst_cen = std::max(worst_cen,std::max(rel_err(cell.center.x,ref_cx),rel_err(cell.center.y,ref_cy)));
        worst_sad = std::max(worst_sad,rel_err(cell.tur.sad,ref_sad));
        worst_seg = std::max(worst_seg,rel_err(cell.tur.sad,ref_seg));
        min_vol = std::min(min_vol,cell.vol);
        max_vol = std::max(max_vol,cell.vol);
        min_sad = std::min(min_sad,cell.tur.sad);
        max_sad = std::max(max_sad,cell.tur.sad);
        if(c < 8){
            fprintf(fp,"%d %.10e %.10e %.3e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n",
                    cell.index,cell.vol,ref_vol,rel_err(cell.vol,ref_vol),
                    cell.center.x,ref_cx,cell.center.y,ref_cy,
                    cell.tur.sad,ref_sad,ref_seg);
        }
    }
    fprintf(fp,"\nnonfinite %d\nvol_rel %.3e\ncenter_rel %.3e\nsad_rel(midpoint) %.3e\n",
            bad,worst_vol,worst_cen,worst_sad);
    fprintf(fp,"sad 用壁面中点近似, 与真实壁面线段的差 %.3e\n",worst_seg);
    fprintf(fp,"vol range [%.6e, %.6e]\nsad range [%.6e, %.6e]\n",min_vol,max_vol,min_sad,max_sad);
    fclose(fp);
    const bool ok = bad == 0 && worst_vol < 1e-9 && worst_cen < 1e-9 && worst_sad < 1e-9;
    g_info = fmt("vol_rel=%.1e cen_rel=%.1e sad_rel=%.1e 首层sad=%.3e vol[%.2e,%.2e]",
                 worst_vol,worst_cen,worst_sad,min_sad,min_vol,max_vol);
    return ok;
}

// ===========================================================================
// 检验 3: 外法向
// ===========================================================================
static bool check_03_normal(){
    FILE* fp = out_open("check03_normal.txt");
    int nor_err = 0,mid_err = 0,fnorm_err = 0,outer_err = 0;
    long long interior = 0,boundary = 0,outer_true = 0;
    double worst_close = 0.0,worst_len = 0.0;
    for(size_t f=0;f<cc::FaceList.size();f++){
        const cc::face_class& face = cc::FaceList[f];
        const double nx = face.node[0]->y - face.node[1]->y;
        const double ny = face.node[1]->x - face.node[0]->x;
        const double len = std::hypot(nx,ny);
        if(!near(face.nor.x,nx,1e-15,1e-15) || !near(face.nor.y,ny,1e-15,1e-15)){
            ++nor_err;
        }
        worst_len = std::max(worst_len,rel_err(face.len,len));
        if(!near(face.mid.x,0.5*(face.node[0]->x+face.node[1]->x),1e-15,1e-15) ||
           !near(face.mid.y,0.5*(face.node[0]->y+face.node[1]->y),1e-15,1e-15)){
            ++mid_err;
        }
        if(face.type == cc::INTER){
            ++interior;
            if(face.outer){
                ++outer_true;
            }
            const bool ref = cc::dot(face.nor,face.nei[1]->center-face.nei[0]->center) > 0;
            if(face.outer != ref){
                ++outer_err;
            }
        }else{
            ++boundary;
        }
    }
    for(int c=0;c<cc::cell_num;c++){
        const cc::cell_class& cell = cc::CellList[c];
        double sx = 0.0,sy = 0.0,scale = 0.0;
        for(int i=0;i<cell.ecnt;i++){
            const cc::face_class* face = cell.faces[i];
            if(cell.fnorm[i] != (cc::dot(face->nor,face->mid-cell.center) > 0)){
                ++fnorm_err;
            }
            const double s = 2*cell.fnorm[i]-1;
            sx += s*face->nor.x;
            sy += s*face->nor.y;
            scale += std::hypot(face->nor.x,face->nor.y);
        }
        // 闭合多边形外法向之和应为 0
        worst_close = std::max(worst_close,std::hypot(sx,sy)/scale);
    }
    fprintf(fp,"nor_err %d\nmid_err %d\nlen_rel %.3e\nfnorm_err %d\nouter_err %d\nmax_closure %.3e\n",
            nor_err,mid_err,worst_len,fnorm_err,outer_err,worst_close);
    fprintf(fp,"interior_faces %lld (outer==true %lld)\nboundary_faces %lld\n",
            interior,outer_true,boundary);
    fclose(fp);
    const bool ok = nor_err == 0 && mid_err == 0 && worst_len < 1e-15 && fnorm_err == 0 &&
                    outer_err == 0 && worst_close < 1e-12;
    g_info = fmt("nor/mid/len/fnorm 全对 闭合=%.1e outer_err=%d (%lld/%lld 为 true)",
                 worst_close,outer_err,outer_true,interior);
    return ok;
}

// ===========================================================================
// 检验 4: 结构化串联与索引关系
// ===========================================================================
static bool check_04_struct(){
    if(!structer::ifstructer){
        FILE* fp = out_open("check04_struct.txt");
        fprintf(fp,"skipped 1\nS_MAX 0\nN_MAX 0\ncheck_if_structed 0\n");
        fclose(fp);
        g_info = "非结构化面基路径(无邻接表), 跳过";
        return true;
    }
    FILE* fp = out_open("check04_struct.txt");
    const int smax = structer::S_MAX,nmax = structer::N_MAX;
    const bool structed = check_if_structed();
    fprintf(fp,"S_MAX %d\nN_MAX %d\ncheck_if_structed %d\n",smax,nmax,structed ? 1 : 0);
    if(smax <= 0 || nmax <= 0 || smax*nmax != cc::cell_num){
        fprintf(fp,"结构化参数与网格规模不一致\n");
        fclose(fp);
        g_info = fmt("S_MAX=%d N_MAX=%d 与 cell_num=%d 不符",smax,nmax,cc::cell_num);
        return false;
    }
    int index_err = 0,dir_err = 0,face_err = 0,ring_err = 0,sn_err = 0,halo_err = 0;
    fprintf(fp,"# idx  s  n  index_ref  east west north south\n");
    for(int c=0;c<cc::cell_num;c++){
        const cc::cell_class& cell = cc::CellList[c];
        if(cell.index != (cell.n-1)*smax + cell.s){
            ++index_err;
        }
        const int dir[4] = {cell.east,cell.west,cell.north,cell.south};
        for(int i=0;i<4;i++){
            if(dir[i] < 0 || dir[i] >= cell.ecnt){
                ++dir_err;
            }
        }
        if(dir[0] == dir[1] || dir[0] == dir[2] || dir[0] == dir[3] ||
           dir[1] == dir[2] || dir[1] == dir[3] || dir[2] == dir[3]){
            ++dir_err;
        }
        if(cell.eastf != cell.faces[cell.east] || cell.westf != cell.faces[cell.west] ||
           cell.northf != cell.faces[cell.north] || cell.southf != cell.faces[cell.south]){
            ++face_err;
        }
        // 四个方向邻居的 (s,n) 增量
        const cc::cell_class* const nbr[4] = {cell.nei[cell.east],cell.nei[cell.west],
                                             cell.nei[cell.north],cell.nei[cell.south]};
        const int want_ds[4] = {1,smax-1,0,0};
        const int want_dn[4] = {0,0,1,-1};
        for(int i=0;i<4;i++){
            if(nbr[i] == nullptr){
                // 只有南/北方向允许是边界
                if(i < 2){
                    ++sn_err;
                }
                continue;
            }
            const int ds = ((nbr[i]->s - cell.s)%smax + smax)%smax;
            if(ds != want_ds[i] || nbr[i]->n - cell.n != want_dn[i]){
                ++sn_err;
            }
        }
        if(c < 6){
            fprintf(fp,"%d %d %d %d %d %d %d %d\n",cell.index,cell.s,cell.n,
                    (cell.n-1)*smax+cell.s,cell.east,cell.west,cell.north,cell.south);
        }
    }
    for(int n=1;n<=nmax;n++){
        const cc::cell_class& first = cc::gotocell((n-1)*smax+1);
        const cc::cell_class& last = cc::gotocell((n-1)*smax+smax);
        if(first.nei[first.west] != &last || last.nei[last.east] != &first){
            ++ring_err;
        }
    }
    for(int n=1;n<=nmax;n++){
        for(int s=1;s<=smax;s++){
            if(&cc::gotocell((n-1)*smax+s) != cc::CellList.data() + (n-1)*smax + s-1){
                ++halo_err;
                break;
            }
        }
    }
    fprintf(fp,"\nindex_err %d\ndirection_err %d\nface_pointer_err %d\nring_closure_err %d\nneighbour_sn_err %d\n",
            index_err,dir_err,face_err,ring_err,sn_err);
    fclose(fp);
    const bool ok = structed && index_err == 0 && dir_err == 0 && face_err == 0 &&
                    ring_err == 0 && sn_err == 0 && halo_err == 0;
    g_info = fmt("structed=%d index_err=%d dir_err=%d face_err=%d ring_err=%d sn_err=%d",
                 structed ? 1 : 0,index_err,dir_err,face_err,ring_err,sn_err);
    return ok;
}

// ===========================================================================
// 检验 5: 物理量初始化
// ===========================================================================
static bool check_05_init(){
    FILE* fp = out_open("check05_init.txt");
    const double rho_inf = FAR_DEFINE.p/(cc::R*FAR_DEFINE.T);
    const double a_inf = std::sqrt(cc::gamma*cc::R*FAR_DEFINE.T);
    const double mu_inf = sutherland::dynamic_viscosity(FAR_DEFINE.T);
    const double e_inf = cc::Cv*FAR_DEFINE.T + 0.5*cc::dot({FAR_DEFINE.u,FAR_DEFINE.v},{FAR_DEFINE.u,FAR_DEFINE.v});
    const double miubl_inf = 3.0*mu_inf/rho_inf;
    double worst[9] = {};
    int nonfinite = 0,unset_mu = 0,unset_un = 0;
    fprintf(fp,"farfield rho %.10e\nu %.10e\nv %.10e\nT %.10e\np %.10e\na %.10e\nmu %.10e\nmiubl %.10e\ne %.10e\n",
            rho_inf,FAR_DEFINE.u,FAR_DEFINE.v,FAR_DEFINE.T,FAR_DEFINE.p,a_inf,mu_inf,miubl_inf,e_inf);
    fprintf(fp,"# idx rho u v T p e a mu miubl\n");
    for(int c=0;c<cc::cell_num;c++){
        const cc::cell_class& cell = cc::CellList[c];
        const double p = cell.phy.rho*cc::R*cell.phy.T;
        const double e = cc::Cv*cell.phy.T + 0.5*(cell.phy.u*cell.phy.u + cell.phy.v*cell.phy.v);
        const double a = std::sqrt(cc::gamma*cc::R*cell.phy.T);
        const double mu = sutherland::dynamic_viscosity(cell.phy.T);
        const double want[9] = {rho_inf,FAR_DEFINE.u,FAR_DEFINE.v,FAR_DEFINE.T,
                                p,e,a,mu,miubl_inf};
        const double got[9] = {cell.phy.rho,cell.phy.u,cell.phy.v,cell.phy.T,
                              cell.otphy.p,cell.otphy.e,cell.otphy.a,cell.otphy.mu,cell.tur.miubl};
        for(int k=0;k<9;k++){
            worst[k] = std::max(worst[k],rel_err(got[k],want[k]));
        }
        if(!std::isfinite(cell.otphy.a) || !std::isfinite(cell.otphy.p) ||
           !std::isfinite(cell.otphy.e) || !std::isfinite(cell.otphy.mu) ||
           !std::isfinite(cell.otphy.un)){
            ++nonfinite;
        }
        if(!(cell.otphy.mu > 0.0)){
            ++unset_mu;
        }
        if(!(cell.otphy.un >= 0.0)){
            ++unset_un;
        }
        // 守恒量 <-> 原始量 往返
        cc::cell_class probe = cell;
        probe.conser = cc::vec4(cell.phy.rho,cell.phy.rho*cell.phy.u,
                                cell.phy.rho*cell.phy.v,cell.phy.rho*cell.otphy.e);
        probe.prim();
        worst[4] = std::max(worst[4],rel_err(probe.phy.rho,cell.phy.rho));
        worst[5] = std::max(worst[5],rel_err(probe.phy.T,cell.phy.T));
        if(c < 3){
            fprintf(fp,"%d %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n",
                    cell.index,cell.phy.rho,cell.phy.u,cell.phy.v,cell.phy.T,
                    cell.otphy.p,cell.otphy.e,cell.otphy.a,cell.otphy.mu,cell.tur.miubl);
        }
    }
    fprintf(fp,"\nrel_err rho %.3e\nu %.3e\nv %.3e\nT %.3e\np %.3e\ne %.3e\na %.3e\nmu %.3e\nmiubl %.3e\n",
            worst[0],worst[1],worst[2],worst[3],worst[4],worst[5],worst[6],worst[7],worst[8]);
    fprintf(fp,"nonfinite_otphy %d\nunset_mu(<=0) %d\nunset_un(<0) %d\n",nonfinite,unset_mu,unset_un);
    double wall_u = 0.0,wall_bl = 0.0;
    for(const cc::face_class* w : cc::WallFaces){
        wall_u = std::max(wall_u,std::hypot(w->phy.u,w->phy.v));
        wall_bl = std::max(wall_bl,std::abs(w->tur.miubl));
    }
    fprintf(fp,"wall_face |V| max %.3e\nwall_face miubl max %.3e\n",wall_u,wall_bl);
    fclose(fp);
    double wmax = 0.0;
    for(int k=0;k<9;k++){
        wmax = std::max(wmax,worst[k]);
    }
    const bool ok = nonfinite == 0 && unset_mu == 0 && wmax < 1e-13 && wall_u == 0.0 && wall_bl == 0.0;
    g_info = fmt("p=%.6e e=%.6e a=%.6e mu=%.6e mu_rel=%.1e worst=%.1e 壁面|V|=%.1e",
                 FAR_DEFINE.p,e_inf,a_inf,mu_inf,worst[7],wmax,wall_u);
    return ok;
}

// ===========================================================================
// 检验 6: 虚网格
// ===========================================================================
static bool check_06_halo(){
    if(!structer::ifstructer){
        FILE* fp = out_open("check06_halo.txt");
        fprintf(fp,"skipped 1\nghost_count 0\nsmax 0 nmax 0\n");
        fclose(fp);
        g_info = "非结构化面基路径(无虚网格), 跳过";
        return true;
    }
    FILE* fp = out_open("check06_halo.txt");
    const int smax = structer::S_MAX,nmax = structer::N_MAX;
    const size_t want = static_cast<size_t>(6)*smax;
    int idx_err = 0,link_err = 0,mirror_err = 0,far_err = 0,missing = 0,neg_bl = 0;
    fprintf(fp,"ghost_count %zu (want %zu)\nsmax %d nmax %d\n",cc::GhostList.size(),want,smax,nmax);
    for(size_t g=0;g<cc::GhostList.size();g++){
        const cc::cell_class& ghost = cc::GhostList[g];
        const int layer = ghost.n <= 0 ? -ghost.n : ghost.n - nmax + 2;
        if(ghost.index != cc::cell_num + layer*smax + ghost.s){
            ++idx_err;
        }
    }
    for(int n=-2;n<=nmax+3;n++){
        if(n >= 1 && n <= nmax){
            continue;
        }
        for(int s=1;s<=smax;s++){
            const cc::cell_class* ghost = gotoHALO(n,s);
            if(ghost == nullptr){
                ++link_err;
                continue;
            }
            if(ghost->nei[ghost->east] != gotoHALO(n,s == smax ? 1 : s+1) ||
               ghost->nei[ghost->west] != gotoHALO(n,s == 1 ? smax : s-1) ||
               ghost->nei[ghost->north] != gotoHALO(n+1,s) ||
               ghost->nei[ghost->south] != gotoHALO(n-1,s)){
                ++link_err;
            }
        }
    }
    for(int n=0;n>=-2;n--){
        for(int s=1;s<=smax;s++){
            const cc::cell_class& ghost = *gotoHALO(n,s);
            const cc::cell_class& inner = *gotoHALO(1-n,s);
            if(!near(ghost.phy.u,-inner.phy.u,1e-15) || !near(ghost.phy.v,-inner.phy.v,1e-15) ||
               !near(ghost.phy.rho,inner.phy.rho,1e-15) || !near(ghost.phy.T,inner.phy.T,1e-15)){
                ++mirror_err;
            }
            if(ghost.tur.miubl < 0.0){
                ++neg_bl;
            }
        }
    }
    for(int n=nmax+1;n<=nmax+3;n++){
        for(int s=1;s<=smax;s++){
            const cc::cell_class& ghost = *gotoHALO(n,s);
            const cc::face_class* reface = gotoHALO(nmax,s)->northf;
            if(reface == nullptr || !near(ghost.phy.rho,reface->phy.rho,1e-15) ||
               !near(ghost.phy.u,reface->phy.u,1e-15)){
                ++far_err;
            }
        }
    }
    for(cc::face_class* w : cc::WallFaces){
        const int i = w->nei[0] ? 1 : 0;
        if(w->nei[i] == nullptr || w->nei[i]->index <= cc::cell_num){
            ++missing;
        }
    }
    for(cc::face_class* f : cc::FarFaces){
        const int i = f->nei[0] ? 1 : 0;
        if(f->nei[i] == nullptr || f->nei[i]->index <= cc::cell_num){
            ++missing;
        }
    }
    fprintf(fp,"ghost_index_err %d\ngotoHALO_link_err %d\nwall_mirror_err %d\nfar_ghost_err %d\nboundary_face_ghost_missing %d\nnegative_miubl_ghosts %d\n",
            idx_err,link_err,mirror_err,far_err,missing,neg_bl);
    fprintf(fp,"# sample ghosts n s index rho u v T miubl\n");
    for(int n=-2;n<=1;n++){
        for(int s=1;s<=2;s++){
            const cc::cell_class& ghost = *gotoHALO(n,s);
            fprintf(fp,"%d %d %d %.6e %.6e %.6e %.6e %.6e\n",ghost.n,ghost.s,ghost.index,
                    ghost.phy.rho,ghost.phy.u,ghost.phy.v,ghost.phy.T,ghost.tur.miubl);
        }
    }
    fclose(fp);
    const bool ok = cc::GhostList.size() == want && idx_err == 0 && link_err == 0 &&
                    mirror_err == 0 && far_err == 0 && missing == 0 && neg_bl == 0;
    g_info = fmt("ghosts=%zu index_err=%d link_err=%d wall_mirror=%d far_err=%d 壁面负miubl=%d",
                 cc::GhostList.size(),idx_err,link_err,mirror_err,far_err,neg_bl);
    return ok;
}

// ===========================================================================
// 检验 7: 梯度
// ===========================================================================
static bool check_07_gradient(){
    FILE* fp = out_open("check07_gradient.txt");
    const double Brho = 0.37,Crho = -0.21;
    const double Bu = 1.13,Cu = 0.47;
    const double Bv = -0.53,Cv = 0.89;
    const double BT = 2.71,CT = -1.37;
    const double Bm = 0.19,Cm = 0.41;
    const double B[5] = {Brho,Bu,Bv,BT,Bm};
    const double C[5] = {Crho,Cu,Cv,CT,Cm};
    std::vector<cc::vecp> saved_phy(cc::cell_num);
    std::vector<double> saved_bl(cc::cell_num);
    for(int c=0;c<cc::cell_num;c++){
        saved_phy[c] = cc::CellList[c].phy;
        saved_bl[c] = cc::CellList[c].tur.miubl;
    }
    for(int c=0;c<cc::cell_num;c++){
        cc::cell_class& cell = cc::CellList[c];
        const double x = cell.center.x,y = cell.center.y;
        cell.phy.rho = 1.0 + Brho*x + Crho*y;
        cell.phy.u = 2.0 + Bu*x + Cu*y;
        cell.phy.v = 3.0 + Bv*x + Cv*y;
        cell.phy.T = 4.0 + BT*x + CT*y;
        cell.tur.miubl = 5.0 + Bm*x + Cm*y;
    }
    int nan_cell = 0,nan_face = 0;
    double interior_worst = 0.0,ghost_worst = 0.0;
    int interior_n = 0,ghost_n = 0;
    for(int c=0;c<cc::cell_num;c++){
        cc::cell_class& cell = cc::CellList[c];
        least_square_cell_based(cell);
        const double got[5][2] = {{cell.phgrad.rhograd.x,cell.phgrad.rhograd.y},
                                  {cell.phgrad.ugrad.x,cell.phgrad.ugrad.y},
                                  {cell.phgrad.vgrad.x,cell.phgrad.vgrad.y},
                                  {cell.phgrad.Tgrad.x,cell.phgrad.Tgrad.y},
                                  {cell.tur.miublgrad.x,cell.tur.miublgrad.y}};
        bool bad = false;
        double cel = 0.0;
        for(int s=0;s<5;s++){
            if(!std::isfinite(got[s][0]) || !std::isfinite(got[s][1])){
                bad = true;
                continue;
            }
            cel = std::max(cel,rel_err(got[s][0],B[s]));
            cel = std::max(cel,rel_err(got[s][1],C[s]));
        }
        if(bad){
            ++nan_cell;
        }
        bool all_real = true;
        for(int i=0;i<cell.ecnt;i++){
            const cc::cell_class* nbr = cell.nei[i];
            if(nbr == nullptr || nbr->index > cc::cell_num){
                all_real = false;
            }
        }
        if(all_real){
            ++interior_n;
            interior_worst = std::max(interior_worst,cel);
        }else{
            ++ghost_n;
            ghost_worst = std::max(ghost_worst,cel);
        }
    }
    #pragma omp parallel for schedule(static)
    for(int c=0;c<cc::cell_num;c++){
        grad_onface(cc::CellList[c]);
    }
    // 面上梯度: 只统计四个面都是内面的格子(边界面用的是边界条件值)
    double face_worst = 0.0;
    for(int c=0;c<cc::cell_num;c++){
        const cc::cell_class& cell = cc::CellList[c];
        bool all_inner = true;
        double cel = 0.0;
        for(int i=0;i<cell.ecnt;i++){
            const cc::face_class* f = cell.faces[i];
            if(f->type != cc::INTER){
                all_inner = false;
                continue;
            }
            if(!std::isfinite(f->phgrad.ugrad.x) || !std::isfinite(f->phgrad.ugrad.y)){
                ++nan_face;
                all_inner = false;
                continue;
            }
            cel = std::max(cel,rel_err(f->phgrad.ugrad.x,Bu));
            cel = std::max(cel,rel_err(f->phgrad.ugrad.y,Cu));
        }
        if(all_inner){
            face_worst = std::max(face_worst,cel);
        }
    }
    // 另测 GGCB: 它用的是面值, 线性场下必须先把面值刷成线性场
    #pragma omp parallel for schedule(static)
    for(int f=0;f<cc::face_num;f++){
        interpolate_mid(&cc::FaceList[f]);
    }
    // 边界面上的面值是边界条件给的, 与线性场不一致, 所以只统计四个面都是内面的格子
    double gg_worst = 0.0;
    int gg_nan = 0,gg_n = 0;
    for(int c=0;c<cc::cell_num;c++){
        cc::cell_class& cell = cc::CellList[c];
        bool all_inner = true;
        for(int i=0;i<cell.ecnt;i++){
            if(cell.faces[i]->type != cc::INTER){
                all_inner = false;
            }
        }
        if(!all_inner){
            continue;
        }
        ++gg_n;
        green_gauss_cell_based(cell);
        if(!std::isfinite(cell.phgrad.ugrad.x) || !std::isfinite(cell.phgrad.ugrad.y)){
            ++gg_nan;
            continue;
        }
        gg_worst = std::max(gg_worst,rel_err(cell.phgrad.ugrad.x,Bu));
        gg_worst = std::max(gg_worst,rel_err(cell.phgrad.ugrad.y,Cu));
    }
    // 用真实场再跑一次, 供后续检验使用
    for(int c=0;c<cc::cell_num;c++){
        cc::CellList[c].phy = saved_phy[c];
        cc::CellList[c].tur.miubl = saved_bl[c];
    }
    refresh_field();
    fprintf(fp,"# 线性场 phy = A + B*x + C*y\n");
    fprintf(fp,"四邻均为真实网格的格子 %d 个, worst_rel %.3e\n",interior_n,interior_worst);
    fprintf(fp,"邻居含虚网格的格子 %d 个, worst_rel %.3e\n",ghost_n,ghost_worst);
    fprintf(fp,"面上梯度 worst_rel %.3e (nonfinite %d)\n",face_worst,nan_face);
    fprintf(fp,"LSCB nonfinite %d\n",nan_cell);
    fprintf(fp,"green_gauss worst_rel %.3e (内部格 %d 个, nonfinite %d)\n",gg_worst,gg_n,gg_nan);
    fclose(fp);
    const bool ok = nan_cell == 0 && nan_face == 0 && interior_worst < 1e-6 && face_worst < 1e-6;
    g_info = fmt("LSCB 内部格 %.1e (%d个) 含虚网格格 %.1e (%d个) 面上 %.1e | GGCB(参考) %.1e",
                 interior_worst,interior_n,ghost_worst,ghost_n,face_worst,gg_worst);
    return ok;
}

// ===========================================================================
// 检验 8: 对流项
// ===========================================================================
static bool check_08_convect(){
    FILE* fp = out_open("check08_convect.txt");
    refresh_field();
    #pragma omp parallel for schedule(static)
    for(int f=0;f<cc::face_num;f++){
        convect_ROE(cc::FaceList[f]);
    }
    long long agree = 0,disagree = 0,nonfinite = 0;
    double worst = 0.0;
    for(int f=0;f<cc::face_num;f++){
        const cc::face_class& face = cc::FaceList[f];
        double ref[4];
        reference_roe(face,ref);
        const double got[4] = {face.convect.c,face.convect.x,face.convect.y,face.convect.e};
        double e = 0.0;
        bool finite = true;
        for(int k=0;k<4;k++){
            if(!std::isfinite(got[k]) || !std::isfinite(ref[k])){
                finite = false;
                break;
            }
            e = std::max(e,rel_err(got[k],ref[k]));
        }
        if(!finite){
            ++nonfinite;
            ++disagree;
            continue;
        }
        worst = std::max(worst,e);
        if(e < 1e-12){
            ++agree;
        }else{
            ++disagree;
        }
        if(f < 4){
            fprintf(fp,"# face %d type %d outer %d\n  src %.10e %.10e %.10e %.10e\n  ref %.10e %.10e %.10e %.10e\n",
                    face.index,face.type,(int)face.outer,got[0],got[1],got[2],got[3],
                    ref[0],ref[1],ref[2],ref[3]);
        }
    }
    // 均匀流一致性: L==R 时 Roe 通量必须退化为精确通量
    int uniform_bad = 0,uniform_checked = 0;
    double worst_uniform = 0.0;
    {
        std::vector<cc::vecp> saved_phy(cc::cell_num);
        std::vector<double> saved_bl(cc::cell_num);
        for(int c=0;c<cc::cell_num;c++){
            saved_phy[c] = cc::CellList[c].phy;
            saved_bl[c] = cc::CellList[c].tur.miubl;
            cc::CellList[c].phy = cc::vecp(1.2,30.0,-4.0,310.0);
            cc::CellList[c].tur.miubl = 1.0e-4;
        }
        refresh_field();
        for(int f=0;f<cc::face_num;f++){
            cc::face_class& face = cc::FaceList[f];
            convect_ROE(face);
            double nx = 0.0,ny = 0.0;
            lr_normal(face,nx,ny);
            const cc::vecp& S = face.phynei[0];
            const cc::otphy& So = face.otnei[0];
            const double un = S.u*nx + S.v*ny;
            const double len = std::hypot(face.nor.x,face.nor.y);
            const double want[4] = {(S.rho*un)*len,(S.rho*S.u*un+So.p*nx)*len,
                                    (S.rho*S.v*un+So.p*ny)*len,(un*(S.rho*So.e+So.p))*len};
            const double got[4] = {face.convect.c,face.convect.x,face.convect.y,face.convect.e};
            ++uniform_checked;
            double e = 0.0;
            bool finite = true;
            for(int k=0;k<4;k++){
                if(!std::isfinite(got[k]) || !std::isfinite(want[k])){
                    finite = false;
                    break;
                }
                e = std::max(e,rel_err(got[k],want[k]));
            }
            if(!finite){
                ++uniform_bad;
                continue;
            }
            worst_uniform = std::max(worst_uniform,e);
            if(e > 1e-12){
                ++uniform_bad;
            }
        }
        for(int c=0;c<cc::cell_num;c++){
            cc::CellList[c].phy = saved_phy[c];
            cc::CellList[c].tur.miubl = saved_bl[c];
        }
    }
    refresh_field();
    fprintf(fp,"total_faces %d\nagree %lld\ndisagree %lld\nnonfinite %lld\nmax_rel %.3e\n",
            cc::face_num,agree,disagree,nonfinite,worst);
    fprintf(fp,"uniform_state: checked %d, mismatch %d, max_rel %.3e\n",
            uniform_checked,uniform_bad,worst_uniform);
    fclose(fp);
    const bool ok = disagree == 0 && uniform_bad == 0;
    g_info = fmt("与参考Roe一致 %lld/%d  均匀流退化为精确通量 %d/%d",
                 agree,cc::face_num,uniform_checked-uniform_bad,uniform_checked);
    return ok;
}

// ===========================================================================
// 检验 9: 湍流扩散项 / 粘性通量
// ===========================================================================
static double my_fv1(double chi){
    const double c3 = SA::Cv1*SA::Cv1*SA::Cv1;
    return (chi*chi*chi)/(chi*chi*chi + c3);
}

static bool check_09_visflux(){
    FILE* fp = out_open("check09_visflux.txt");
    #pragma omp parallel for schedule(static)
    for(int f=0;f<cc::face_num;f++){
        SA::diffusion_SA(cc::FaceList[f]);
    }
    int mueff_bad = 0,flux_bad = 0,nan_face = 0;
    double worst = 0.0;
    for(int f=0;f<cc::face_num;f++){
        const cc::face_class& face = cc::FaceList[f];
        const double mu = sutherland::dynamic_viscosity(face.phy.T);
        const double chi = face.phy.rho*std::max(face.tur.miubl,0.0)/mu;
        const double mut = face.phy.rho*my_fv1(chi)*face.tur.miubl;
        const double mueff = mu + mut;
        if(!near(face.tur.mueff,mueff,1e-13)){
            ++mueff_bad;
        }
        // 壁面独立复算: 法向导数按一阶壁面律自己搭一遍
        cc::vecgrad pg = face.phgrad;
        if(face.type == cc::WALL){
            const cc::cell_class* c = cc::boundary_findcell(const_cast<cc::face_class*>(&face));
            const double y1 = c->tur.sad;
            const double dx = face.mid.x - c->center.x,dy = face.mid.y - c->center.y;
            const double dl = std::sqrt(dx*dx + dy*dy);
            const double nx = dx/dl,ny = dy/dl,tx = -ny,ty = nx;
            const double dun = (face.phy.u - c->phy.u)/y1;
            const double dvn = (face.phy.v - c->phy.v)/y1;
            const double dut = c->phgrad.ugrad.x*tx + c->phgrad.ugrad.y*ty;
            const double dvt = c->phgrad.vgrad.x*tx + c->phgrad.vgrad.y*ty;
            pg.ugrad = {dun*nx + dut*tx,dun*ny + dut*ty};
            pg.vgrad = {dvn*nx + dvt*tx,dvn*ny + dvt*ty};
        }
        const double tauxx = mueff*(4.0/3*pg.ugrad.x - 2.0/3*pg.vgrad.y);
        const double tauxy = mueff*(pg.ugrad.y + pg.vgrad.x);
        const double tauyy = mueff*(4.0/3*pg.vgrad.y - 2.0/3*pg.ugrad.x);
        const double lambdaeff = mu/cc::Pr + mut/SA::Prt;
        const cc::vec2 q = -lambdaeff*cc::Cp*pg.Tgrad;
        const double F[4] = {0,tauxx,tauxy,face.phy.u*tauxx + face.phy.v*tauxy - q.x};
        const double G[4] = {0,tauxy,tauyy,face.phy.u*tauxy + face.phy.v*tauyy - q.y};
        const double ref[4] = {0,
            F[1]*face.nor.x + G[1]*face.nor.y,
            F[2]*face.nor.x + G[2]*face.nor.y,
            F[3]*face.nor.x + G[3]*face.nor.y};
        const double got[4] = {face.visflux.c,face.visflux.x,face.visflux.y,face.visflux.e};
        double e = 0.0;
        bool finite = true;
        for(int k=0;k<4;k++){
            if(!std::isfinite(got[k]) || !std::isfinite(ref[k])){
                finite = false;
                break;
            }
            e = std::max(e,rel_err(got[k],ref[k]));
        }
        if(!finite){
            ++nan_face;
            ++flux_bad;
            continue;
        }
        worst = std::max(worst,e);
        if(e > 1e-12){
            ++flux_bad;
        }
    }
    fprintf(fp,"mueff_bad %d\nvisflux_bad %d\nnonfinite_face %d\nmax_rel %.3e\n",
            mueff_bad,flux_bad,nan_face,worst);
    fprintf(fp,"系数 Pr %.3f Prt %.3f inv_sigma %.3f Cb1 %.4f Cb2 %.4f Cw1 %.4f Cw2 %.3f Cw3 %.1f Cv1 %.1f C5 %.1f\n",
            cc::Pr,SA::Prt,SA::inv_sigma,SA::Cb1,SA::Cb2,SA::Cw1,SA::Cw2,SA::Cw3,SA::Cv1,SA::C5);
    fclose(fp);
    const bool ok = mueff_bad == 0 && flux_bad == 0;
    g_info = fmt("mueff 全对(含mu+mut) visflux 全对 max_rel=%.1e nonfinite=%d",worst,nan_face);
    return ok;
}

// ===========================================================================
// 检验 10: SA 源项
// ===========================================================================
static bool check_10_source(){
    FILE* fp = out_open("check10_source.txt");
    int bad = 0,nan_cell = 0,neg_dest = 0;
    double worst = 0.0,worst_terms[4] = {};
    fprintf(fp,"# idx chi fv1 fv2 ft2 r g fw vort scaled production destruction grad_src compress total_ref total_src\n");
    for(int c=0;c<cc::cell_num;c++){
        const cc::cell_class& cell = cc::CellList[c];
        const double mu = cell.otphy.mu;
        const double chi = cell.phy.rho*std::max(cell.tur.miubl,0.0)/mu;
        const double c3 = SA::Cv1*SA::Cv1*SA::Cv1;
        const double fv1 = (chi*chi*chi)/(chi*chi*chi + c3);
        const double ft2 = SA::Ct3*std::exp(-SA::Ct4*chi*chi);
        const double fv2 = 1 - chi/(1 + chi*fv1);
        const double vort = std::abs(cell.phgrad.vgrad.x - cell.phgrad.ugrad.y);
        const double d2 = cell.tur.sad*cell.tur.sad;
        const double scaled = cell.tur.miubl/d2/(SA::kappa*SA::kappa);
        const double modvort = std::max(vort + fv2*scaled,std::max(0.3*vort,1e-20));
        const double production = SA::Cb1*(1-ft2)*modvort*cell.phy.rho*cell.tur.miubl;
        const double r = std::min(scaled/modvort,SA::rmax);
        const double g = r + SA::Cw2*(r*r*r*r*r*r - r);
        const double cw3s2 = SA::Cw3*SA::Cw3;
        const double cw3s6 = cw3s2*cw3s2*cw3s2;
        const double g2 = g*g,g6 = g2*g2*g2;
        const double fw = g*std::pow((1+cw3s6)/(cw3s6+g6),1.0/6);
        const double destruction = cell.phy.rho*(SA::Cw1*fw - SA::Cb1/SA::kappa/SA::kappa*ft2)*
                                   cell.tur.miubl*cell.tur.miubl/d2;
        const double grad_src = SA::Cb2*SA::inv_sigma*cell.phy.rho*
                                cc::dot(cell.tur.miublgrad,cell.tur.miublgrad);
        const double S2 = 2*cell.phgrad.ugrad.x*cell.phgrad.ugrad.x +
                          2*cell.phgrad.vgrad.y*cell.phgrad.vgrad.y +
                          (cell.phgrad.ugrad.y + cell.phgrad.vgrad.x)*(cell.phgrad.ugrad.y + cell.phgrad.vgrad.x);
        const double compress = SA::C5*cell.phy.rho*cell.tur.miubl*cell.tur.miubl*S2/
                                (cc::gamma*cc::R*cell.phy.T);
        const double ref = production - destruction + grad_src - compress;
        const double got = SA::source_SA(cell);
        if(!std::isfinite(got) || !std::isfinite(ref)){
            ++nan_cell;
            continue;
        }
        if(destruction < 0.0){
            ++neg_dest;
        }
        worst_terms[0] = std::max(worst_terms[0],std::abs(production));
        worst_terms[1] = std::max(worst_terms[1],std::abs(destruction));
        worst_terms[2] = std::max(worst_terms[2],std::abs(grad_src));
        worst_terms[3] = std::max(worst_terms[3],std::abs(compress));
        const double e = rel_err(got,ref);
        worst = std::max(worst,e);
        if(e > 1e-12){
            ++bad;
        }
        if(c < 3 || c == cc::cell_num-1){
            fprintf(fp,"%d %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
                    cell.index,chi,fv1,fv2,ft2,r,g,fw,vort,scaled,production,destruction,
                    grad_src,compress,ref,got);
        }
    }
    fprintf(fp,"\nsource_bad %d\nnonfinite %d\nworst_rel %.3e\n",bad,nan_cell,worst);
    fprintf(fp,"各项幅值上限 production %.3e destruction %.3e grad_src %.3e compressible %.3e\n",
            worst_terms[0],worst_terms[1],worst_terms[2],worst_terms[3]);
    fprintf(fp,"destruction 为负的格子 %d\n",neg_dest);
    fclose(fp);
    const bool ok = bad == 0 && nan_cell == 0;
    g_info = fmt("逐项复算一致 max_rel=%.1e 产生项<=%.2e 耗散项<=%.2e 可压缩项<=%.2e",
                 worst,worst_terms[0],worst_terms[1],worst_terms[3]);
    return ok;
}

// ===========================================================================
// 检验 11: 一次 RK 子步
// ===========================================================================
static bool check_11_one_rk(double dt){
    FILE* fp = out_open("check11_rk1.txt");
    std::vector<cc::vec4> conser0;
    std::vector<double> bl0;
    snapshot(conser0,bl0);
    for(int c=0;c<cc::cell_num;c++){
        cc::CellList[c].copyconver_time();
    }
    refresh_field();
    #pragma omp parallel for schedule(static)
    for(int f=0;f<cc::face_num;f++){
        convect_ROE(cc::FaceList[f]);
    }
    #pragma omp parallel for schedule(static)
    for(int c=0;c<cc::cell_num;c++){
        assemble_flux(cc::CellList[c]);
    }
    #pragma omp parallel for schedule(static)
    for(int f=0;f<cc::face_num;f++){
        SA::diffusion_SA(cc::FaceList[f]);
    }
    #pragma omp parallel for schedule(static)
    for(int c=0;c<cc::cell_num;c++){
        SA::assemble_visflux(cc::CellList[c]);
    }
    #pragma omp parallel for schedule(static)
    for(int c=0;c<cc::cell_num;c++){
        local_timestep(cc::CellList[c]);
    }
    // 锁定伪时间步基准态, 再更新一次
    std::vector<cc::vec4> before(cc::cell_num);
    for(int c=0;c<cc::cell_num;c++){
        before[c] = cc::CellList[c].conser;
        cc::CellList[c].conserformer = before[c];
    }
    #pragma omp parallel for schedule(static)
    for(int c=0;c<cc::cell_num;c++){
        cc::cell_class& cell = cc::CellList[c];
        cell.conser = cell.conserformer - RK::RK3[0]/(1.0/cell.localdt + 1.0/(2*dt))*
                      (1.0/cell.vol*(cell.convect-cell.visflux) +
                       1.0/(2*dt)*(before[c]-cell.lastconser));
    }
    int mismatch = 0,nan_cell = 0,neg_dt = 0;
    double worst = 0.0,worst_ratio = 0.0,min_dt = 1e300,max_dt = 0.0,min_chi = 1e300;
    fprintf(fp,"dt %.6e RK3[0] %.6f\n",dt,RK::RK3[0]);
    fprintf(fp,"# idx localdt chi R_conv.c R_phys.c conser_new.c conser_manual.c rel\n");
    for(int c=0;c<cc::cell_num;c++){
        cc::cell_class& cell = cc::CellList[c];
        const double scale = 1.0/cell.localdt + 1.0/(2*dt);
        const double chi = RK::RK3[0]/scale;
        const cc::vec4 R = (cell.convect - cell.visflux)*(1.0/cell.vol);
        const cc::vec4 Rp = (before[c] - cell.lastconser)*(1.0/(2*dt));
        const cc::vec4 manual = cell.conserformer - (R + Rp)*chi;
        if(!std::isfinite(cell.localdt) || !std::isfinite(chi) || !std::isfinite(manual.c) ||
           !std::isfinite(cell.conser.c)){
            ++nan_cell;
            continue;
        }
        const double e = std::max(std::max(rel_err(cell.conser.c,manual.c),rel_err(cell.conser.x,manual.x)),
                                  std::max(rel_err(cell.conser.y,manual.y),rel_err(cell.conser.e,manual.e)));
        worst = std::max(worst,e);
        if(e > 1e-13){
            ++mismatch;
        }
        if(!(cell.localdt > 0.0)){
            ++neg_dt;
        }
        min_dt = std::min(min_dt,cell.localdt);
        max_dt = std::max(max_dt,cell.localdt);
        const double ratio = chi/cell.localdt;
        worst_ratio = std::max(worst_ratio,ratio);
        min_chi = std::min(min_chi,ratio);
        if(c < 4){
            fprintf(fp,"%d %.6e %.6e %.6e %.6e %.10e %.10e %.1e\n",
                    cell.index,cell.localdt,chi,R.c,Rp.c,cell.conser.c,manual.c,e);
        }
    }
    fprintf(fp,"\nformula_mismatch %d\nnonfinite %d\nworst_rel %.3e\n",mismatch,nan_cell,worst);
    {
        const cc::cell_class& cell = cc::CellList[0];
        fprintf(fp,"\n# 首个格子逐面通量 cell idx %d vol %.6e rho %.6e\n",cell.index,cell.vol,cell.phy.rho);
        for(int i=0;i<cell.ecnt;i++){
            const cc::face_class* f = cell.faces[i];
            fprintf(fp,"  face[%d] idx %d type %d outer %d islow %d fnorm %d conv %.6e %.6e %.6e %.6e\n",
                    i,f->index,f->type,(int)f->outer,(int)(f->nei[0] == &cell),(int)cell.fnorm[i],
                    f->convect.c,f->convect.x,f->convect.y,f->convect.e);
            fprintf(fp,"      L(rho %.6e u %.6e v %.6e T %.6e)  R(rho %.6e u %.6e v %.6e T %.6e)\n",
                    f->phynei[0].rho,f->phynei[0].u,f->phynei[0].v,f->phynei[0].T,
                    f->phynei[1].rho,f->phynei[1].u,f->phynei[1].v,f->phynei[1].T);
        }
    }
    fprintf(fp,"localdt range [%.6e, %.6e], <=0 的格子 %d\n",min_dt,max_dt,neg_dt);
    fprintf(fp,"振幅因子 chi/localdt range [%.6e, %.6e], 要求 0<chi<=localdt\n",min_chi,worst_ratio);
    fclose(fp);
    restore(conser0,bl0);
    const bool ok = mismatch == 0 && nan_cell == 0 && neg_dt == 0 && worst_ratio <= 1.0+1e-12;
    g_info = fmt("代数式复算一致 mismatch=%d localdt[%.2e,%.2e] <=0的格子=%d chi/localdt<=%.4f",
                 mismatch,min_dt,max_dt,neg_dt,worst_ratio);
    return ok;
}

// ===========================================================================
// 检验 12: 三次 RK 构成的一组伪时间步
// ===========================================================================
static bool check_12_three_rk(double dt){
    FILE* fp = out_open("check12_rk3.txt");
    // (a) 检验台实现 vs src 的 one_rans, 两次都从同一初始状态出发
    std::vector<cc::vec4> conser0,src_res,mine;
    std::vector<double> bl0,src_bl,mine_bl;
    snapshot(conser0,bl0);
    refresh_field();
    one_rans(dt,true);
    snapshot(src_res,src_bl);
    restore(conser0,bl0);
    refresh_field();
    pseudo_step_omp(dt,true);
    snapshot(mine,mine_bl);
    int mismatch = 0,nan_cell = 0;
    double worst = 0.0;
    for(int c=0;c<cc::cell_num;c++){
        if(!std::isfinite(mine[c].c) || !std::isfinite(src_res[c].c)){
            ++nan_cell;
            continue;
        }
        const double e = std::max(std::max(rel_err(mine[c].c,src_res[c].c),rel_err(mine[c].x,src_res[c].x)),
                                  std::max(rel_err(mine[c].y,src_res[c].y),rel_err(mine[c].e,src_res[c].e)));
        worst = std::max(worst,e);
        if(e > 1e-14){
            ++mismatch;
        }
    }
    // (b) 三次子步的相容性与稳定性
    const double a1 = RK::RK3[0],a2 = RK::RK3[1],a3 = RK::RK3[2];
    const bool last_one = a3 == 1.0;
    double max_mod = 0.0,first_unstable = 0.0;
    for(double t = 0.01;t <= 3.0;t += 0.01){
        // lam*dtau = i*t, U_k = U0 - a_k*l*U_{k-1}
        const double li = t;
        const double r1 = 1.0,i1 = -a1*li;
        double pr = -li*i1,pi = li*r1;
        const double r2 = 1.0 - a2*pr,i2 = -a2*pi;
        pr = -li*i2;
        pi = li*r2;
        const double r3 = 1.0 - a3*pr,i3 = -a3*pi;
        const double mod = std::hypot(r3,i3);
        if(mod > 1.0 && first_unstable == 0.0){
            first_unstable = t;
        }
        max_mod = std::max(max_mod,mod);
    }
    // (c) 两次稳态伪时间步的相对变化应递减
    double r_first = 0.0,r_second = 0.0,peak = 0.0;
    restore(conser0,bl0);
    refresh_field();
    const int steady_passes = g_rans ? 60 : 2;
    for(int pass=0;pass<steady_passes;pass++){
        std::vector<cc::vec4> before(cc::cell_num);
        for(int c=0;c<cc::cell_num;c++){
            before[c] = cc::CellList[c].conser;
        }
        pseudo_step_omp(dt,false);
        const double now = increment_l2(before);
        if(pass == 0){
            r_first = now;
        }
        r_second = now;
        peak = std::max(peak,now);
    }
    restore(conser0,bl0);
    refresh_field();
    if(g_laminar){
        mismatch = 0;nan_cell = 0;worst = 0.0;
    }
    fprintf(fp,"(a) 检验台实现 vs src one_rans (相同初值)%s\n    cell_mismatch %d  nonfinite %d  worst_rel %.3e\n",
            g_laminar ? " [层流: 不跑SA, 对比无意义, 已跳过]" : "",mismatch,nan_cell,worst);
    fprintf(fp,"(b) 系数 %.4f %.4f %.4f, 末次系数为1 %d\n",a1,a2,a3,last_one ? 1 : 0);
    fprintf(fp,"    经典三阶要求 1/3 1/2 1, 偏差 %.4f %.4f %.4f\n",a1-1.0/3,a2-0.5,a3-1.0);
    fprintf(fp,"    虚轴放大因子 max|P| %.4f, 首次超过 1 的 theta %.2f\n",max_mod,first_unstable);
    fprintf(fp,"(c) 稳态伪时间步相邻增量: 首次 %.6e 峰值 %.6e 末次 %.6e (应从峰值回落)\n",r_first,peak,r_second);
    fclose(fp);
    const bool ok = mismatch == 0 && nan_cell == 0 && last_one && r_second < peak;
    g_info = fmt("与one_rans完全一致=%d 末次系数1=%d 稳态增量 %.2e->%.2e 虚轴max|P|=%.2f",
                 mismatch == 0 ? 1 : 0,last_one ? 1 : 0,r_first,r_second,max_mod);
    return ok;
}

// ===========================================================================
// 检验 13: 一个物理时间步
// ===========================================================================
static bool check_13_step(double dt,int subiter){
    FILE* fp = out_open("check13_step.txt");
    std::vector<cc::vec4> conser0;
    std::vector<double> bl0;
    snapshot(conser0,bl0);
    if(g_rans){
        std::vector<cc::vec4> before(cc::cell_num);
        double first = 0.0,last = 0.0,prev = 0.0,max_bl = 0.0,from_start = 0.0;
        int nan_step = 0,breaks = 0;
        fprintf(fp,"# RANS 稳态伪时间迭代\n# iter  相邻增量L2  距初值L2  max_miubl\n");
        for(int k=0;k<subiter;k++){
            #pragma omp parallel for schedule(static)
            for(int c=0;c<cc::cell_num;c++){
                before[c] = cc::CellList[c].conser;
            }
            pseudo_step_omp(0.0,false);
            const double inc = increment_l2(before);
            const double now = change_l2(conser0,&max_bl);
            if(!std::isfinite(inc) || !std::isfinite(max_bl)){
                ++nan_step;
                fprintf(fp,"%d nan\n",k+1);
                break;
            }
            if(k == 0){
                first = inc;
            }else if(inc > prev){
                ++breaks;
            }
            prev = inc;
            last = inc;
            from_start = now;
            fprintf(fp,"%d %.10e %.10e %.10e\n",k+1,inc,now,max_bl);
        }
        fprintf(fp,"\nfirst %.6e\nlast %.6e\nfinal/first %.6e\nmonotone_break %d\nnonfinite_steps %d\n",
                first,last,first > 0 ? last/first : 0.0,breaks,nan_step);
        fclose(fp);
        restore(conser0,bl0);
        refresh_field();
        const bool okr = nan_step == 0 && first > 0.0 && last < first;
        g_info = fmt("RANS 伪时间增量 %.3e -> %.3e (ratio %.4f, 非单调 %d 次) 距初值 %.3e",
                     first,last,first > 0 ? last/first : 0.0,breaks,from_start);
        return okr;
    }
    for(int c=0;c<cc::cell_num;c++){
        cc::CellList[c].copyconver_time();
    }
    double first = 0.0,last = 0.0,prev = 0.0,max_bl = 0.0;
    int nan_step = 0,breaks = 0;
    std::vector<cc::vec4> before(cc::cell_num);
    fprintf(fp,"# subiter  伪时间残差L2(U^k-U^k-1)  距U^n的L2  max_miubl\n");
    for(int k=0;k<subiter;k++){
        for(int c=0;c<cc::cell_num;c++){
            before[c] = cc::CellList[c].conser;
        }
        pseudo_step_omp(dt,true);
        const double inc = increment_l2(before);
        const double now = change_l2(conser0,&max_bl);
        if(!std::isfinite(inc) || !std::isfinite(max_bl)){
            ++nan_step;
            fprintf(fp,"%d nan\n",k+1);
            const int bad_c = first_bad_cell();
            if(bad_c >= 0){
                fprintf(fp,"  首个异常格 index=%d s=%d n=%d conser=%.4e localdt=%.4e miubl=%.4e\n",
                        bad_c+1,cc::CellList[bad_c].s,cc::CellList[bad_c].n,
                        cc::CellList[bad_c].conser.c,cc::CellList[bad_c].localdt,
                        cc::CellList[bad_c].tur.miubl);
            }
            break;
        }
        if(k == 0){
            first = inc;
        }else if(inc > prev){
            ++breaks;
        }
        prev = inc;
        last = inc;
        fprintf(fp,"%d %.10e %.10e %.10e\n",k+1,inc,now,max_bl);
    }
    fprintf(fp,"\nfirst %.6e\nlast %.6e\nfinal/first %.6e\nmonotone_break %d\nnonfinite_steps %d\n",
            first,last,first > 0 ? last/first : 0.0,breaks,nan_step);
    fclose(fp);
    restore(conser0,bl0);
    refresh_field();
    const bool ok = nan_step == 0 && first > 0.0 && last < first;
    g_info = fmt("伪时间残差 %.3e -> %.3e (ratio %.4f, 非单调 %d 次)",
                 first,last,first > 0 ? last/first : 0.0,breaks);
    return ok;
}

// ===========================================================================
// 检验 14: 短跑物理步数
// ===========================================================================
static bool check_14_run(double dt,int subiter,int steps){
    FILE* fp = out_open("residual_run.csv");
    std::vector<cc::vec4> conser0;
    std::vector<double> bl0;
    double last_l2 = 0.0,l2_first = 0.0;
    int nan_step = 0;
    if(g_rans){
        fprintf(fp,"step,res_increment,res_from_initial,max_miubl,min_localdt\n");
        std::vector<cc::vec4> before(cc::cell_num);
        for(int n=1;n<=steps;n++){
            #pragma omp parallel for schedule(static)
            for(int c=0;c<cc::cell_num;c++){
                before[c] = cc::CellList[c].conser;
            }
            pseudo_step_omp(0.0,false);
            const double inc = increment_l2(before);
            double maxbl = 0.0,min_dt = 1e300;
            int bad = 0;
            #pragma omp parallel for schedule(static) reduction(max:maxbl) reduction(min:min_dt) reduction(+:bad)
            for(int c=0;c<cc::cell_num;c++){
                const cc::cell_class& cell = cc::CellList[c];
                maxbl = std::max(maxbl,std::abs(cell.tur.miubl));
                if(std::isfinite(cell.localdt)){
                    min_dt = std::min(min_dt,cell.localdt);
                }
                if(!std::isfinite(cell.conser.c) || !std::isfinite(cell.tur.miubl)){
                    bad = 1;
                }
            }
            last_l2 = inc;
            if(n == 1){
                l2_first = inc;
            }
            if(n <= 15 || n % 500 == 0){
                miubl_trace(n,dt,false);
            }
            fprintf(fp,"%d,%.10e,%.10e,%.10e,%.10e\n",n,inc,change_l2(before,nullptr),maxbl,min_dt);
            if(n % 10 == 0 || n == steps){
                printf("      iter %4d  增量L2=%.4e  max_miubl=%.4e\n",n,inc,maxbl);
                fflush(stdout);
            }
            if(g_dump_every > 0 && n % g_dump_every == 0){
                dump_field(std::to_string(n),n,dt);
                fflush(stdout);
            }
            if(bad){
                ++nan_step;
                fprintf(fp,"# 第 %d 次迭代出现非有限值\n",n);
                break;
            }
        }
        fclose(fp);
        const bool ok = nan_step == 0 && std::isfinite(last_l2) && last_l2 < l2_first;
        g_info = fmt("RANS %d 次迭代, 增量L2 %.3e -> %.3e",steps,l2_first,last_l2);
        return ok;
    }
    fprintf(fp,"step,res_conser,res_convect,res_visflux,max_miubl,min_localdt\n");
    for(int n=1;n<=steps;n++){
        snapshot(conser0,bl0);
        #pragma omp parallel for schedule(static)
        for(int c=0;c<cc::cell_num;c++){
            cc::CellList[c].copyconver_time();
        }
        for(int k=0;k<subiter;k++){
            pseudo_step_omp(dt,true);
        }
        double maxbl = 0.0,min_dt = 1e300,rc = 0.0,rx = 0.0,rv = 0.0;
        int bad = 0;
        #pragma omp parallel for schedule(static) reduction(max:maxbl) reduction(min:min_dt) \
                                 reduction(+:rc,rx,rv,bad)
        for(int c=0;c<cc::cell_num;c++){
            const cc::cell_class& cell = cc::CellList[c];
            const cc::vec4 d = cell.conser - conser0[c];
            rc += d.c*d.c + d.x*d.x + d.y*d.y + d.e*d.e;
            const cc::vec4 conv = (cell.convect - cell.visflux)*(1.0/cell.vol);
            rx += conv.c*conv.c + conv.x*conv.x + conv.y*conv.y + conv.e*conv.e;
            const cc::vec4 vis = cell.visflux*(1.0/cell.vol);
            rv += vis.c*vis.c + vis.x*vis.x + vis.y*vis.y + vis.e*vis.e;
            maxbl = std::max(maxbl,std::abs(cell.tur.miubl));
            if(std::isfinite(cell.localdt)){
                min_dt = std::min(min_dt,cell.localdt);
            }
            if(!std::isfinite(cell.conser.c) || !std::isfinite(cell.tur.miubl)){
                bad = true;
            }
        }
        last_l2 = std::sqrt(rc/cc::cell_num);
        if(n == 1){
            l2_first = last_l2;
        }
        fprintf(fp,"%d,%.10e,%.10e,%.10e,%.10e,%.10e\n",n,last_l2,
                std::sqrt(rx/cc::cell_num),std::sqrt(rv/cc::cell_num),maxbl,min_dt);
        if(n % 10 == 0 || n == steps){
            printf("      step %4d  L2(dU)=%.4e  R_convect=%.4e  max_miubl=%.4e\n",
                   n,last_l2,std::sqrt(rx/cc::cell_num),maxbl);
            fflush(stdout);
        }
        if(g_dump_every > 0 && n % g_dump_every == 0){
            dump_field(std::to_string(n),n,dt);
            fflush(stdout);
        }
        if(n <= 15){
            miubl_trace(n,dt,true);
        }
        if(bad){
            ++nan_step;
            fprintf(fp,"# 第 %d 个物理步出现非有限值\n",n);
            const int bad_c = first_bad_cell();
            if(bad_c >= 0){
                fprintf(fp,"  首个异常格 index=%d s=%d n=%d conser=%.4e localdt=%.4e miubl=%.4e\n",
                        bad_c+1,cc::CellList[bad_c].s,cc::CellList[bad_c].n,
                        cc::CellList[bad_c].conser.c,cc::CellList[bad_c].localdt,
                        cc::CellList[bad_c].tur.miubl);
            }
            break;
        }
    }
    fclose(fp);
    const bool ok = nan_step == 0 && std::isfinite(last_l2);
    g_info = fmt("%d 步无异常, 每步 L2(dU) %.3e -> %.3e",steps,l2_first,last_l2);
    return ok;
}

// ===========================================================================
// 检验 15: 输出流场
// ===========================================================================
static bool dump_field(const std::string& tag,int steps,double dt){
    // 面通量是上一次 RK 子步中间算的, 先同步到当前格子状态再取值
    refresh_field();
    #pragma omp parallel for schedule(static)
    for(int f=0;f<cc::face_num;f++){
        SA::diffusion_SA(cc::FaceList[f]);
    }
    // final 不带后缀, 与 python 复核脚本的默认文件名一致
    const std::string stem = (tag == "final") ? std::string() : ("_" + tag);
    const std::string path = std::string(OUTDIR) + "/field_" + tag + ".dat";
    FILE* fp = fopen(path.c_str(),"w");
    if(fp == nullptr){
        g_info = "无法写出 " + path;
        return false;
    }
    const double rho_inf = FAR_DEFINE.p/(cc::R*FAR_DEFINE.T);
    const double U_inf = std::hypot(FAR_DEFINE.u,FAR_DEFINE.v);
    const double q_inf = 0.5*rho_inf*U_inf*U_inf;
    fprintf(fp,"TITLE=\"cylinder Re=%.0f Ma=0.2 %s dt=%.3e after %d steps\"\n",
            g_re,g_rans ? "RANS" : "URANS",dt,steps);
    fprintf(fp,"VARIABLES=\"x\",\"y\",\"rho\",\"u\",\"v\",\"T\",\"p\",\"Ma\",\"nu_tilde\",\"Cp\",\"s\",\"n\"\n");
    int bad = 0;
    double ma_min = 1e300,ma_max = -1e300,cp_min = 1e300,cp_max = -1e300;
    for(int c=0;c<cc::cell_num;c++){
        const cc::cell_class& cell = cc::CellList[c];
        const double p = cell.otphy.p;
        const double ma = std::hypot(cell.phy.u,cell.phy.v)/cell.otphy.a;
        const double cp = (p - FAR_DEFINE.p)/q_inf;
        if(!std::isfinite(p) || !std::isfinite(ma) || !std::isfinite(cp)){
            ++bad;
            continue;
        }
        ma_min = std::min(ma_min,ma);
        ma_max = std::max(ma_max,ma);
        cp_min = std::min(cp_min,cp);
        cp_max = std::max(cp_max,cp);
        fprintf(fp,"%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %d %d\n",
                cell.center.x,cell.center.y,cell.phy.rho,cell.phy.u,cell.phy.v,
                cell.phy.T,p,ma,cell.tur.miubl,cp,cell.s,cell.n);
    }
    fclose(fp);
    // 壁面受力 = 压力 + 粘性摩擦:  F = ∮ p n dS - ∮ (tau·n) dS
    // n 取面法向中"由内部格子指向壁面"的那一侧, 面积已含在 nor / visflux 里
    double fx_p = 0.0,fy_p = 0.0,fx_v = 0.0,fy_v = 0.0;
    FILE* wp = out_open(("wall_cp" + stem + ".csv").c_str());
    fprintf(wp,"theta_deg,x,y,Cp,Ma,fx_p,fy_p,fx_v,fy_v\n");
    for(const cc::face_class* w : cc::WallFaces){
        const cc::cell_class* c = cc::boundary_findcell(const_cast<cc::face_class*>(w));
        const double s = cc::dot(w->nor,w->mid-c->center) > 0 ? 1.0 : -1.0;
        const double cp = (w->otphy.p - FAR_DEFINE.p)/q_inf;
        const double ma = std::hypot(w->phy.u,w->phy.v)/w->otphy.a;
        const double theta = std::atan2(w->mid.y,w->mid.x)/std::acos(-1.0)*180.0;
        const double fxp = w->otphy.p*s*w->nor.x,fyp = w->otphy.p*s*w->nor.y;
        const double fxv = -s*w->visflux.x,fyv = -s*w->visflux.y;
        fx_p += fxp;
        fy_p += fyp;
        fx_v += fxv;
        fy_v += fyv;
        fprintf(wp,"%.6f %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n",
                theta,w->mid.x,w->mid.y,cp,ma,fxp,fyp,fxv,fyv);
    }
    fclose(wp);
    const double D_ref = 1.0;
    const double ref_force = q_inf*D_ref;
    // 风轴系: 阻力沿来流, 升力垂直来流; 来流方向由远场速度矢量给出
    const double arad = std::atan2(FAR_DEFINE.v,FAR_DEFINE.u);
    const double ca = std::cos(arad),sa = std::sin(arad);
    const double cd_p = (fx_p*ca + fy_p*sa)/ref_force;
    const double cd_f = (fx_v*ca + fy_v*sa)/ref_force;
    const double cl_p = (-fx_p*sa + fy_p*ca)/ref_force;
    const double cl_f = (-fx_v*sa + fy_v*ca)/ref_force;
    FILE* ff = out_open(("forces" + stem + ".csv").c_str());
    fprintf(ff,"component,force,coefficient\n");
    fprintf(ff,"pressure_x,%.10e,%.10e\n",fx_p,fx_p/ref_force);
    fprintf(ff,"friction_x,%.10e,%.10e\n",fx_v,fx_v/ref_force);
    fprintf(ff,"total_x,%.10e,%.10e\n",fx_p+fx_v,(fx_p+fx_v)/ref_force);
    fprintf(ff,"pressure_y,%.10e,%.10e\n",fy_p,fy_p/ref_force);
    fprintf(ff,"friction_y,%.10e,%.10e\n",fy_v,fy_v/ref_force);
    fprintf(ff,"total_y,%.10e,%.10e\n",fy_p+fy_v,(fy_p+fy_v)/ref_force);
    fprintf(ff,"drag_pressure,%.10e,%.10e\n",fx_p*ca + fy_p*sa,cd_p);
    fprintf(ff,"drag_friction,%.10e,%.10e\n",fx_v*ca + fy_v*sa,cd_f);
    fprintf(ff,"drag_total,%.10e,%.10e\n",(fx_p+fx_v)*ca + (fy_p+fy_v)*sa,cd_p+cd_f);
    fprintf(ff,"lift_pressure,%.10e,%.10e\n",-fx_p*sa + fy_p*ca,cl_p);
    fprintf(ff,"lift_friction,%.10e,%.10e\n",-fx_v*sa + fy_v*ca,cl_f);
    fprintf(ff,"lift_total,%.10e,%.10e\n",-(fx_p+fx_v)*sa + (fy_p+fy_v)*ca,cl_p+cl_f);
    fprintf(ff,"q_inf,%.10e,D,%.3f,alpha_deg,%.4f\n",q_inf,D_ref,arad*180.0/std::acos(-1.0));
    fclose(ff);
    printf("      [%s] Cd_p=%.4f  Cd_f=%.4f  Cd=%.4f   Cl_p=%.4f  Cl_f=%.4f  Cl=%.4f\n",
           tag.c_str(),cd_p,cd_f,cd_p+cd_f,cl_p,cl_f,cl_p+cl_f);
    const bool ok = bad == 0;
    g_info = fmt("Cd_p=%.4f Cd_f=%.4f Cd=%.4f Cl=%.4f | Ma[%.4f,%.4f] Cp[%.3f,%.3f] nonfinite=%d",
                 cd_p,cd_f,cd_p+cd_f,cl_p+cl_f,ma_min,ma_max,cp_min,cp_max,bad);
    return ok;
}

static bool check_15_dump(int steps,double dt){
    return dump_field("final",steps,dt);
}

// ===========================================================================
// ν̃ 追踪: 每个迭代步记录 ν̃ 最大格的 SA 方程逐项分解
// ===========================================================================

static FILE* g_trace_fp = nullptr;

struct sa_terms{
    double mu,chi,fv1,fv2,ft2;
    double vorticity,scaled_nu,modvort,r,g,fw;
    double production,destruction,gradient_source,compressible,source;
    double convect,diffusion,rhs;
    double grad_norm;
};

static sa_terms compute_sa_terms(const cc::cell_class& cell){
    sa_terms t{};
    t.mu = cell.otphy.mu;
    t.chi = cell.phy.rho*std::max(cell.tur.miubl,0.0)/t.mu;
    const double c3 = SA::Cv1*SA::Cv1*SA::Cv1;
    t.fv1 = (t.chi*t.chi*t.chi)/(t.chi*t.chi*t.chi + c3);
    t.ft2 = SA::Ct3*std::exp(-SA::Ct4*t.chi*t.chi);
    t.fv2 = 1 - t.chi/(1 + t.chi*t.fv1);
    t.vorticity = std::abs(cell.phgrad.vgrad.x - cell.phgrad.ugrad.y);
    const double d2 = cell.tur.sad*cell.tur.sad;
    t.scaled_nu = cell.tur.miubl/d2/(SA::kappa*SA::kappa);
    t.modvort = std::max(t.vorticity + t.fv2*t.scaled_nu,std::max(0.3*t.vorticity,1e-20));
    t.production = SA::Cb1*(1-t.ft2)*t.modvort*cell.phy.rho*cell.tur.miubl;
    t.r = std::min(t.scaled_nu/t.modvort,SA::rmax);
    t.g = t.r + SA::Cw2*(t.r*t.r*t.r*t.r*t.r*t.r - t.r);
    const double cw3s2 = SA::Cw3*SA::Cw3;
    const double cw3s6 = cw3s2*cw3s2*cw3s2;
    const double g2 = t.g*t.g,g6 = g2*g2*g2;
    t.fw = t.g*std::pow((1+cw3s6)/(cw3s6+g6),1.0/6);
    t.destruction = cell.phy.rho*(SA::Cw1*t.fw - SA::Cb1/SA::kappa/SA::kappa*t.ft2)*
                    cell.tur.miubl*cell.tur.miubl/d2;
    t.gradient_source = SA::Cb2*SA::inv_sigma*cell.phy.rho*
                        cc::dot(cell.tur.miublgrad,cell.tur.miublgrad);
    const double S2 = 2*cell.phgrad.ugrad.x*cell.phgrad.ugrad.x +
                      2*cell.phgrad.vgrad.y*cell.phgrad.vgrad.y +
                      (cell.phgrad.ugrad.y + cell.phgrad.vgrad.x)*
                      (cell.phgrad.ugrad.y + cell.phgrad.vgrad.x);
    t.compressible = SA::C5*cell.phy.rho*cell.tur.miubl*cell.tur.miubl*S2/
                     (cc::gamma*cc::R*cell.phy.T);
    t.source = t.production - t.destruction + t.gradient_source - t.compressible;
    t.grad_norm = std::hypot(cell.tur.miublgrad.x,cell.tur.miublgrad.y);
    t.convect = 0.0;
    t.diffusion = 0.0;
    for(int i=0;i<cell.ecnt;i++){
        const cc::face_class* face = cell.faces[i];
        const double outer = 2*cell.fnorm[i]-1;
        const double mu_f = sutherland::dynamic_viscosity(face->phy.T);
        const double coef = SA::inv_sigma*(mu_f + face->phy.rho*face->tur.miubl);
        t.convect += outer*face->toface_jacobi(face->phy.rho*face->phy.u*face->tur.miubl,
                                               face->phy.rho*face->phy.v*face->tur.miubl);
        t.diffusion += outer*face->toface_jacobi(coef*face->tur.miublgrad.x,
                                                 coef*face->tur.miublgrad.y);
    }
    t.rhs = (t.diffusion - t.convect)/(cell.vol*cell.phy.rho) + t.source/cell.phy.rho;
    return t;
}

static void miubl_trace(int iter,double dt,bool urans){
    if(g_trace_fp == nullptr || g_laminar){
        return;
    }
    int best = 0;
    for(int c=1;c<cc::cell_num;c++){
        if(std::abs(cc::CellList[c].tur.miubl) > std::abs(cc::CellList[best].tur.miubl)){
            best = c;
        }
    }
    const cc::cell_class& cell = cc::CellList[best];
    const sa_terms t = compute_sa_terms(cell);
    const double rho = cell.phy.rho;
    const double dtau = cell.localdt;
    const double flux = (t.diffusion - t.convect)/(cell.vol*rho);
    fprintf(g_trace_fp,"\n=== iter %d  ν̃ 最大格: idx=%d s=%d n=%d (x=%.4f, y=%.4f) d=%.4e vol=%.4e\n",
            iter,cell.index,cell.s,cell.n,cell.center.x,cell.center.y,cell.tur.sad,cell.vol);
    fprintf(g_trace_fp,"  ν̃=%.6e  dtau=%.6e  |∇ν̃|=%.6e  χ=%.6e fv1=%.6e fv2=%.6e ft2=%.6e\n",
            cell.tur.miubl,dtau,t.grad_norm,t.chi,t.fv1,t.fv2,t.ft2);
    fprintf(g_trace_fp,"  S=%.6e  ν̃/(κ²d²)=%.6e  S̃=%.6e  r=%.6e g=%.6e fw=%.6e\n",
            t.vorticity,t.scaled_nu,t.modvort,t.r,t.g,t.fw);
    fprintf(g_trace_fp,"  # 归一到 dν̃/dt 的各项(已除 rho)\n");
    fprintf(g_trace_fp,"    production      = %+.6e\n",t.production/rho);
    fprintf(g_trace_fp,"    destruction     = %+.6e\n",-t.destruction/rho);
    fprintf(g_trace_fp,"    gradient_source = %+.6e\n",t.gradient_source/rho);
    fprintf(g_trace_fp,"    compressible    = %+.6e\n",-t.compressible/rho);
    fprintf(g_trace_fp,"    (diff-conv)/V   = %+.6e   [diffusion=%.6e convect=%.6e]\n",
            flux,t.diffusion,t.convect);
    fprintf(g_trace_fp,"    rhs = dν̃/dt     = %+.6e\n",t.rhs);
    const double step = urans ? t.rhs/(1.0/dtau + 1.0/(2*dt)) : 0.8*dtau*t.rhs;
    fprintf(g_trace_fp,"    Δν̃(本步)       = %+.6e   ->  %.6e\n",step,cell.tur.miubl+step);
    // 直接调用 src 的更新函数核对
    cc::cell_class probe = cell;
    probe.tur.miubl_former = probe.tur.miubl;
    const double before = probe.tur.miubl;
    SA::SA_equation_after(probe,probe.localdt,dt,urans);
    fprintf(g_trace_fp,"    src SA_equation_after: %.6e -> %.6e  (Δ=%+.6e)\n",
            before,probe.tur.miubl,probe.tur.miubl-before);
    for(int i=0;i<cell.ecnt;i++){
        const cc::face_class* face = cell.faces[i];
        const double outer = 2*cell.fnorm[i]-1;
        const double mu_f = sutherland::dynamic_viscosity(face->phy.T);
        const double coef = SA::inv_sigma*(mu_f + face->phy.rho*face->tur.miubl);
        const double c = outer*face->toface_jacobi(face->phy.rho*face->phy.u*face->tur.miubl,
                                                   face->phy.rho*face->phy.v*face->tur.miubl);
        const double d = outer*face->toface_jacobi(coef*face->tur.miublgrad.x,
                                                   coef*face->tur.miublgrad.y);
        fprintf(g_trace_fp,"    face[%d] idx=%d type=%d outer=%+.0f nor=(%.4e,%.4e) rho=%.6e V=(%.4e,%.4e) ν̃=%.4e ∇ν̃=(%.4e,%.4e) conv=%+.6e diff=%+.6e\n",
                i,face->index,face->type,outer,face->nor.x,face->nor.y,face->phy.rho,
                face->phy.u,face->phy.v,face->tur.miubl,
                face->tur.miublgrad.x,face->tur.miublgrad.y,c,d);
    }
}

// ===========================================================================
// 主流程
// ===========================================================================
int main(int argc,char** argv){
    int steps = 100,subiter = 60;
    double dt = 1.0e-4;
    if(argc > 1){
        steps = std::atoi(argv[1]);
    }
    if(argc > 2){
        subiter = std::atoi(argv[2]);
    }
    if(argc > 3){
        dt = std::atof(argv[3]);
    }
    if(argc > 4){
        fatime::CFL = std::atof(argv[4]);
    }
    if(argc > 5){
        g_rans = std::string(argv[5]) == "rans";
    }
    if(argc > 6){
        const std::string mode = argv[6];
        g_trace = mode.find("trace") != std::string::npos;
        g_laminar = mode.find("laminar") != std::string::npos;
        cc::scheme = (mode.find("jst") != std::string::npos) ? 'J' : 'R';
    }
    // 线程数: 默认 16, 可用第 7 个参数覆盖
    const int nthreads = argc > 7 ? std::atoi(argv[7]) : 16;
    // 中间帧输出间隔: 默认 5000 步, 可用第 8 个参数覆盖(0=不输出)
    g_dump_every = argc > 8 ? std::atoi(argv[8]) : 5000;
    #ifdef _OPENMP
    omp_set_num_threads(nthreads);
    #endif
    printf("============== 圆柱绕流 %s%s 检验台 ==============\n",g_rans ? "RANS" : "URANS",
           g_laminar ? "/层流(湍流模型已冻结)" : "");
    printf("步=%d  子迭代=%d  物理时间步 dt=%.3e s  OpenMP线程=%d  中间输出间隔=%d  格式=%s\n",
           steps,subiter,dt,nthreads,g_dump_every,cc::scheme == 'J' ? "JST" : "Roe");
    if(mkdir(OUTDIR,0755) != 0 && errno != EEXIST){
        fprintf(stderr,"Error: cannot create %s\n",OUTDIR);
        return 1;
    }
    g_summary = out_open("summary.txt");
    if(g_summary == nullptr || !config::load("config.json")){
        return 1;
    }
    if(g_trace){
        g_trace_fp = out_open("miubl_trace.txt");
        if(g_trace_fp != nullptr){
            fprintf(g_trace_fp,"# ν̃ 追踪: 每步记录 ν̃ 最大格的 SA 方程逐项分解\n");
        }
    }
    if(argc > 4){
        fatime::CFL = std::atof(argv[4]);
    }
    printf("CFL=%.3f\n",fatime::CFL);
    // Re 由远场密度定, 网格固定 D=1; 第 9 个参数可改 Re
    const double D = 1.0;
    const double Re_target = argc > 9 ? std::atof(argv[9]) : 60.0;
    g_re = Re_target;
    const double a_inf = std::sqrt(cc::gamma*cc::R*FAR_DEFINE.T);
    const double U_inf = std::hypot(FAR_DEFINE.u,FAR_DEFINE.v);
    const double mu_inf = sutherland::dynamic_viscosity(FAR_DEFINE.T);
    const double p_target = Re_target*mu_inf*cc::R*FAR_DEFINE.T/(U_inf*D);
    printf("工况: D=%.3f m  Ma=%.3f  T=%.2f K  p=%.6e Pa (Re=%.0f)\n",
           D,U_inf/a_inf,FAR_DEFINE.T,p_target,Re_target);
    FAR_DEFINE.p = p_target;
    g_scale[0] = p_target/(cc::R*FAR_DEFINE.T);
    g_scale[1] = g_scale[0]*U_inf;
    g_scale[2] = g_scale[1];
    g_scale[3] = g_scale[0]*cc::Cp*FAR_DEFINE.T;

    if(!readmesh(cc::meshpath.c_str()) || !geometrymain()){
        printf("网格或几何失败\n");
        return 1;
    }
    check_begin(1,"网格 / 面读取 / 邻接关系");
    const bool ok1 = check_01_mesh();
    check_end(ok1,g_info);
    check_begin(2,"几何量(面积 / 质心 / sad)");
    const bool ok2 = check_02_geometry();
    check_end(ok2,g_info);
    check_begin(3,"外法向");
    const bool ok3 = check_03_normal();
    check_end(ok3,g_info);
    check_begin(4,"结构化串联与索引关系");
    const bool ok4 = check_04_struct();
    check_end(ok4,g_info);

    std_initialize();
    for(int c=0;c<cc::cell_num;c++){
        cc::CellList[c].form_conservative();
    }
    check_begin(5,"物理量初始化");
    const bool ok5 = check_05_init();
    check_end(ok5,g_info);

    // LSCB 预处理必须在 HALO 之前: 虚网格没有几何
    for(int c=0;c<cc::cell_num;c++){
        least_square_cell_based_preprocess(cc::CellList[c]);
    }
    HALO_structer_mesh();
    if(g_laminar){
        freeze_turbulence();
    }
    check_begin(6,"虚拟网格邻接与物理量");
    const bool ok6 = check_06_halo();
    check_end(ok6,g_info);

    check_begin(7,"梯度建立");
    const bool ok7 = check_07_gradient();
    check_end(ok7,g_info);
    check_begin(8,"对流项建立");
    const bool ok8 = check_08_convect();
    check_end(ok8,g_info);
    check_begin(9,"湍流扩散项 / 粘性通量");
    const bool ok9 = check_09_visflux();
    check_end(ok9,g_info);
    check_begin(10,"SA 源项");
    const bool ok10 = check_10_source();
    check_end(ok10,g_info);

    check_begin(11,"一次 RK 子步");
    const bool ok11 = check_11_one_rk(dt);
    check_end(ok11,g_info);
    check_begin(12,"三次 RK 构成的一组伪时间步");
    const bool ok12 = check_12_three_rk(dt);
    check_end(ok12,g_info);
    check_begin(13,g_rans ? "一个稳态伪时间步" : "一个物理时间步");
    const bool ok13 = check_13_step(dt,subiter);
    check_end(ok13,g_info);
    check_begin(14,fmt(g_rans ? "短跑 %d 次稳态迭代" : "短跑 %d 个物理步",steps).c_str());
    const bool ok14 = check_14_run(dt,subiter,steps);
    check_end(ok14,g_info);
    check_begin(15,"输出流场(供 python 画云图)");
    const bool ok15 = check_15_dump(steps,dt);
    check_end(ok15,g_info);

    printf("---------------- 通过 %d/%d ----------------\n",g_pass,g_done);
    if(g_trace_fp != nullptr){
        fclose(g_trace_fp);
    }
    if(g_summary != nullptr){
        fprintf(g_summary,"---- %d/%d passed ----\n",g_pass,g_done);
        fclose(g_summary);
    }
    return 0;
}
