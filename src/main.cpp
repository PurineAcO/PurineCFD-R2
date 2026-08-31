#include <cstdio>
#include "classconfig.h"
#include "readmesh.h"
#include "geometry.h"
#include "initialize.h"
#include "config.h"
#include "boundary.h"
#include "convect.h"
#include "interpolat.h"
#include "grad.h"
#include "udf.hpp"

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
    }

    // 第1处函数钩子:定义边界条件
    DEFINE_BOUNDARY();

    // 初始化
    std_initialize();

    // 求解期
    // 这里还有一个巨大的步长/时间推进的循环
    // for(int i=1;i<=cc::max_step;i++){}
    // 这里还有一个RK循环
    // for(int j=1;j<=cc::RK;j++){}
        wall_boundary();velocity_inlet_boundary();pressure_outlet_boundary();    // 三个边界条件
        for(cc::cell_class& cell : cc::CellList){
            green_gauss_cell_based(cell);   // 建立梯度
            interpolate_mid(cell);          // 做面上的重构(基本量),含边界
            cell.form_physic();               // 建立cell的全部物理量
            cell.form_conservative();         // 建立守恒量
            convect_JST(cell);              // 建立对流项            
        }

    return 0;
}