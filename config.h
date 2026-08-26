#pragma once

#include "classconfig.h"
namespace cc {

    inline int cell_num = 0;                        // 网格总数
    inline int face_num = 0;                        // 面总数
    inline const int HALO = 2;                      // HALO网格层数(可能弃用)
    inline const char* path = "mesh/tunnel.txt";    // 网格文件位置

    // 安全的网格访问
    cell_class& gotocell(int number);

    // 安全的面访问
    face_class& gotoface(int number);
}