#include <cstdio>
#include <cmath>
#include "classconfig.h"
#include "readmesh.h"
#include "geometry.h"
#include "initialize.h"
#include "config.h"
#include "boundary.h"
#include "convect.h"
#include "interpolat.h"
#include "dissipation.h"
#include "timarch.h"
#include "residual.h"
#include "SA.h"
#include "io.h"

#define allcell for(cc::cell_class& cell : cc::CellList)

int main(){
    config::load();
    freopen(cc::testpath, "w", stdout);
    if(readmesh(cc::meshpath)){return 1;}
    geometrymain();
    allcell sad(cell);
    std_initialize();
    for(cc::cell_class& cell : cc::CellList){ cell.form_conservative(); }
    dump_field(0);

    // 定常: 当地时间步长迭代至残差收敛
    const int dump_step = config::dump_step;   // 场输出间隔
    const int conv_check = config::conv_step;  // 残差检查间隔
    int step = 0;
    double res_init = -1.0;
    double sim_t = 0.0;
    while(step < cc::max_step && (cc::total_time <= 0.0 || sim_t < cc::total_time)){
        step++;
        allcell cell.copyconver();
        allcell local_timestep(cell);
        if(fatime::USE_GLOBAL_DT){
            double dt = 1e30;
            for(const auto& c : cc::CellList){ if(c.localdt < dt) dt = c.localdt; }
            for(auto& c : cc::CellList){ c.localdt = dt; }
            sim_t += dt;
        }
        for(int j=0;j<5;j++){
            slip_wall_boundary(); far_field_boundary();
            allcell cell.reform();
            allcell interpolate_mid(cell);
            allcell cell.form_physic();
            allcell jst::shockwave_recognize(cell);
            allcell jst::laplace_dissipation(cell);
            allcell cell.form_conservative();
            allcell convect_JST(cell);
            allcell SA::diffusion_SA(cell);
            allcell jst::JST_dissipation(cell);
            for(auto& cell : cc::CellList){
                for(int s=0;s<4;s++){
                    double R = (cell.convect[s] - cell.disspiation.Fd[s] - cell.tur.Ft[s])/cell.vol;
                    cell.conser[s] = cell.conserformer[s] - RK::RK[j]*cell.localdt*R;
                }
            }
            if(cc::viscous){
                for(auto& cell : cc::CellList){ SA::SA_equation_RK(cell, RK::RK[j]); }
            }
        }
        res::report_update(step);
        {
            bool nan_found = false;
            for(const cc::cell_class& c : cc::CellList){
                if(!(std::isfinite(c.phy.rho)&&std::isfinite(c.phy.u)&&
                     std::isfinite(c.phy.v)&&std::isfinite(c.phy.T)) ||
                   !std::isfinite(c.tur.miubl)){
                    printf("[NaN] step %d first bad cell#%d (%.6f,%.6f) rho=%.3e miubl=%.3e\n",
                           step, c.index, c.center.x, c.center.y, c.phy.rho, c.tur.miubl);
                    nan_found = true; break;
                }
            }
            if(nan_found){ printf("NaN detected at step %d, abort.\n", step); break; }
        }
        if(step % dump_step == 0){ dump_field(step); save_checkpoint(step); }
        if(step % conv_check == 0){
            double res = res::residual_absmax();
            if(res_init < 0){ res_init = res; }
            if(res < res_init*1e-4){
                printf("Converged at step %d, res=%.3e\n", step, res);
                break;
            }
        }
    }
    dump_field(step);
    save_checkpoint(step);
    printf("Total step: %d\n", step);
    return 0;
}
