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
    BEFORE_CONFIG();    // 可能会被弃用
    freopen(cc::testpath, "w", stdout);
    if(readmesh(cc::meshpath)){return 1;}
    geometrymain();
    DEFINE_BOUNDARY();  // 可能会被弃用
    std_initialize();
    for(cc::cell_class& cell : cc::CellList){ cell.form_conservative(); }
    dump_field(0);

    // 定常: 当地时间步长迭代至残差收敛
    const int dump_step = 10000;    // 每N步存一次场
    const int conv_check = 200;     // 每N步查残差
    int step = 0;
    double res_init = -1.0;
    while(step < cc::max_step){
        step++;
        allcell cell.copyconver();
        allcell local_timestep(cell);
        for(int j=0;j<5;j++){
            slip_wall_boundary(); far_field_boundary();
            allcell cell.reform();
            allcell interpolate_mid(cell);
            allcell cell.form_physic();
            allcell jst::shockwave_recognize(cell);
            allcell jst::laplace_dissipation(cell);
            allcell cell.form_conservative();
            allcell convect_JST(cell);
            allcell jst::JST_dissipation(cell);
            for(auto& cell : cc::CellList){
                for(int s=0;s<4;s++){
                    double R = (cell.convect[s] - cell.disspiation.Fd[s])/cell.vol;
                    cell.conser[s] = cell.conserformer[s] - RK::RK[j]*cell.localdt*R;
                }
            }
        }
        res::report_update(step);
        if(step % dump_step == 0){ dump_field(step); save_checkpoint(step); }
        if(step % conv_check == 0){
            double res = res::current_residual_all();
            if(res_init < 0){ res_init = res; }
            if(res < res_init*1e-4){       // 残差下降4个量级
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
