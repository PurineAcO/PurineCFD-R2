#include <cstdio>
#include <cmath>
#include <vector>
#include "classconfig.h"
#include "readmesh.h"
#include "geometry.h"
#include "initialize.h"
#include "config.h"
#include "boundary.h"
#include "convect.h"
#include "interpolat.h"
#include "grad.h"
#include "SA.h"
#include "dissipation.h"
#include "timarch.h"
#include "residual.h"
#include "udf.hpp"

#define allcell for(cc::cell_class& cell : cc::CellList)

int main(){

    // 第0处函数钩子:定义宏
    BEFORE_CONFIG();

    // 测试输出文件路径
    freopen(cc::testpath,"w",stdout);

    // 网格读取和生成期
    if(readmesh(cc::meshpath)){return 1;}
    if(linkmesh()){return 1;}

    // 网格几何的生成期
    for(cc::cell_class& cell: cc::CellList){
        findnode(cell);
        volume(cell);
        center(cell);
        walldistance_alternative(cell,cc::H);
        // printf("%f\n",cell.sad);
        cell.face_normal_out();
        jst::JST_dissipation_INIT(cell);
    }

    // 第1处函数钩子:定义边界条件
    DEFINE_BOUNDARY();

    // 初始化
    std_initialize();
    wall_boundary();velocity_inlet_boundary();pressure_outlet_boundary();
    for(cc::cell_class& cell: cc::CellList){cell.form_conservative();}

    // 记录初始场作为第0步
    res::report_update(0);

    // 求解期
    for(int step=1;step<=cc::max_step;step++){
    allcell cell.copyconver();
    allcell local_timestep(cell);
    if(fatime::USE_GLOBAL_DT){
        double dtmin = 1e30;
        for(auto& c : cc::CellList){ if(c.localdt < dtmin){ dtmin = c.localdt; } }
        for(auto& c : cc::CellList){ c.localdt = dtmin; }
    }
    // RK循环
    for(int j=0;j<5;j++){
        wall_boundary();velocity_inlet_boundary();pressure_outlet_boundary();
        allcell cell.reform();
        allcell green_gauss_cell_based(cell);
        allcell interpolate_mid(cell);
        allcell cell.form_physic();
        allcell jst::shockwave_recognize(cell);
        allcell jst::laplace_dissipation(cell);
        allcell cell.form_conservative();
        allcell convect_JST(cell);
        allcell SA::SA_diffusion(cell);
        allcell SA::SA_source(cell);
        allcell jst::JST_dissipation(cell);
        allcell {std::vector<double> conser_RK;conser_RK.resize(5);
                for(int s=0;s<5;s++){
                conser_RK[s] = (cell.convect[s] - cell.diffusion[s] - cell.disspiation.Fd[s])/cell.vol - cell.source[s];
                cell.conser[s] = cell.conserformer[s] - RK::RK[j] * cell.localdt * conser_RK[s];}
            }
        }
        res::report_update(step);
    }
    return 0;
}