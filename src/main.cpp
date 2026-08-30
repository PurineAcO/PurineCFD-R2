#include <cstdio>
#include "classconfig.h"
#include "readmesh.h"
#include "geometry.h"
#include "initialize.h"
#include "config.h"
#include "udf.hpp"

int main(){

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

    // 第一处函数钩子:定义边界条件
    DEFINE_BOUNDARY();

    // 初始化
    std_initialize();

    return 0;
}