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
#include "io.h"
#include "udf.hpp"

#define allcell for(cc::cell_class& cell : cc::CellList)

int main(){
    BEFORE_CONFIG();
    // freopen(cc::testpath, "w", stdout);

    if(readmesh(cc::meshpath)){ return 1; }
    if(linkmesh()){ return 1; }

    for(cc::cell_class& cell : cc::CellList){
        findnode(cell);
        volume(cell);
        center(cell);
        cell.face_normal_out();
    }

    DEFINE_BOUNDARY();

    double ckpt_t = 0.0;
    bool resumed = load_checkpoint(ckpt_t);
    if(resumed){
        cc::total_time = ckpt_t;
        printf("From time: %.6f go on\n", cc::total_time);
    }else{
        cc::total_time = 0.0;
        std_initialize();
    }
    slip_wall_boundary(); velocity_inlet_boundary(); pressure_outlet_boundary();
    for(cc::cell_class& cell : cc::CellList){ cell.form_conservative(); }
    if(!resumed){ dump_field(0.0); }

    res::report_update(0);

    const double t_end  = cc::total_time + 0.10;
    const double dt_dump = 0.004;
    double next_dump = std::floor(cc::total_time/dt_dump)*dt_dump + dt_dump;
    int step = 0;
    while(cc::total_time < t_end && step < cc::max_step){
        step++;
        allcell cell.copyconver();
        allcell local_timestep(cell);
        double dtmin = 1e30;
        for(auto& c : cc::CellList){ if(c.localdt < dtmin){ dtmin = c.localdt; } }
        for(auto& c : cc::CellList){ c.localdt = dtmin; }
        for(int j=0;j<5;j++){
            slip_wall_boundary(); velocity_inlet_boundary(); pressure_outlet_boundary();
            allcell cell.reform();
            allcell interpolate_mid(cell);
            allcell cell.form_physic();
            allcell jst::shockwave_recognize(cell);
            allcell jst::laplace_dissipation(cell);
            allcell cell.form_conservative();
            allcell convect_JST(cell);
            allcell jst::JST_dissipation(cell);
            for(auto& cell : cc::CellList){
                for(int s=0;s<cc::NEQ;s++){
                    double rk = (cell.convect[s] - cell.disspiation.Fd[s])/cell.vol;
                    cell.conser[s] = cell.conserformer[s] - RK::RK[j]*cell.localdt*rk;
                }
            }
        }
        cc::total_time += dtmin;
        res::report_update(step);
        if(cc::total_time >= next_dump){
            dump_field(cc::total_time);
            save_checkpoint(cc::total_time);
            next_dump += dt_dump;
        }
    }
    save_checkpoint(cc::total_time);
    printf("Total step: %d ,total time:%.6f s\n", step, cc::total_time);
    return 0;
}
