#pragma once

#include "classconfig.h"
namespace cc {

    // 求解器设置
    inline int cell_num = 0;                        // 网格总数
    inline int face_num = 0;                        // 面总数
    inline int type_total = 0;                      // 边界条件类型总数
    inline const int HALO = 2;                      // HALO网格层数(可能弃用)
    inline const char* path = "mesh/tunnel.txt";    // 网格文件位置

    // 网格和边界条件约定
    inline const short INTER = 0;                     // 内部面
    inline const short WALL = 1;                      // 无滑移壁面
    inline const short POL = 2;                       // 压力出口
    inline const short VIL = 3;                       // 速度入口
    inline const short FAR = 4;                       // 压力远场  

    // 安全的网格访问
    cell_class& gotocell(int number);

    // 安全的面访问
    face_class& gotoface(int number);

    // 判断面的类型
    short get_facetype(const char* face_std_name);
}