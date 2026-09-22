#include <cmath>
#include <cstdio>
#include <omp.h>
#include "config.hpp"
#include "classconfig.hpp"
#include "parallel.hpp"
#include "readmesh.hpp"
#include "HALO.hpp"
#include "geometry.hpp"
#include "initialize.hpp"
#include "boundary.hpp"
#include "interpolate.hpp"
#include "grad.hpp"
#include "convect.hpp"
#include "dissipation.hpp"
#include "SA.hpp"
#include "timarch.hpp"
#include "residual.hpp"
#include "io.hpp"

#define allcell for(cc::cell_class& cell : cc::CellList)

static bool check_field(const char* func,int step){
    if(cc::field_bad_cell.load(std::memory_order_relaxed) == 0){
        return true;
    }
    for(const cc::cell_class& cell : cc::CellList){
        if(cc::field_ok(cell)){
            continue;
        }
        fprintf(stderr,
                "Error: Invalid flow state at step %d, function %s, cell #%d (%.6f,%.6f)"
                " rho=%.6e u=%.6e v=%.6e T=%.6e p=%.6e miubl=%.6e\n",
                step,func,cell.index,cell.center.x,cell.center.y,
                cell.phy.rho,cell.phy.u,cell.phy.v,cell.phy.T,cell.phy.p,cell.tur.miubl);
        return false;
    }
    return true;
}

static void recover_flow(){
#pragma omp parallel for schedule(static)
    for(int i=0;i<cc::cell_num;i++){
        cc::CellList[i].reform();
        cc::CellList[i].form_physic();
    }
}

static bool rk_stage(double rk,int step){
    recover_flow();
    if(!check_field("recover_flow",step)){
        return false;
    }
    update_ghost_field();
    slip_wall_boundary();
    far_field_boundary();
#pragma omp parallel for schedule(static)
    for(int i=0;i<cc::cell_num;i++){
        interpolate_mid(cc::CellList[i]);
    }
#pragma omp parallel for schedule(static)
    for(int i=0;i<cc::cell_num;i++){
        cc::cell_class& cell = cc::CellList[i];
        green_gauss_cell_based(cell);
        jst::shockwave_recognize(cell);
        jst::laplace_dissipation(cell);
    }
#pragma omp parallel for schedule(static)
    for(int i=0;i<cc::face_num;i++){
        cc::face_class& face = cc::FaceList[i];
        face.form_physic();
        face_gradient(face);
        convect_JST(face);
        SA::diffusion_SA(face);
    }
#pragma omp parallel for schedule(static)
    for(int i=0;i<cc::cell_num;i++){
        cc::cell_class& cell = cc::CellList[i];
        assemble_flux(cell);
        jst::JST_dissipation(cell);
        SA::SA_equation_RK(cell,rk);
    }
#pragma omp parallel for schedule(static)
    for(int i=0;i<cc::cell_num;i++){
        cc::cell_class& cell = cc::CellList[i];
        for(int s=0;s<4;s++){
            double R = (cell.convect[s] - cell.diss.Fd[s] - cell.visflux[s])*cell.invvol;
            cell.conser[s] = cell.conserformer[s] - rk*cell.localdt*R;
        }
        cell.tur.miubl = cell.tur.miubl_next;
        field_mark(cell);
    }
    if(!check_field("rk_stage",step)){
        return false;
    }
    return true;
}

static bool solve(){
    int step = 0;
    int dumped = 0;
    bool converged = false;
    if(!dump_field(0)){
        return false;
    }
    while(step < cc::max_step){
        step++;
#pragma omp parallel for schedule(static)
        for(int i=0;i<cc::cell_num;i++){
            cc::CellList[i].copyconver();
            local_timestep(cc::CellList[i]);
        }
        for(int j=0;j<5;j++){
            if(!rk_stage(RK::RK[j],step)){
                return false;
            }
        }
        recover_flow();
        if(!check_field("recover_flow",step)){
            return false;
        }
        res::report_update(step);
        if(step % config::dump_step == 0){
            if(!dump_field(step)){
                return false;
            }
            dumped = step;
        }
        if(step % config::conv_step == 0){
            double drop = res::worst_drop();
            if(drop <= res::drop_target){
                printf("Converged at step %d, every field dropped at least %.1f decades\n",
                       step,-std::log10(drop));
                converged = true;
                break;
            }
        }
    }
    if(dumped != step && !dump_field(step)){
        return false;
    }
    if(!converged){
        printf("Iteration limit reached; convergence criterion not satisfied.\n");
    }
    printf("Total step: %d\n",step);
    return true;
}

int main(int argc,char** argv){
    if(argc > 2){
        fprintf(stderr,"Error: Usage: purinecfd [config.json]\n");
        return 1;
    }
    const char* thread_policy = parallel::configure_threads();
    if(!config::load(argc == 2 ? argv[1] : "config.json")){
        return 1;
    }
    if(!open_log(cc::testpath.c_str())){
        return 1;
    }
    printf("Steady SA-RANS | CFL=%.2f | max_steps=%lld | OpenMP threads=%d | thread_policy=%s\n",
           fatime::CFL,cc::max_step,omp_get_max_threads(),thread_policy);
    if(!readmesh(cc::meshpath.c_str())){
        return 1;
    }
    if(!geometrymain()){
        return 1;
    }
    allcell sad(cell);
    std_initialize();
    allcell cell.form_conservative();
    HALO_structer_mesh();
    if(!solve()){
        return 1;
    }
    return 0;
}
