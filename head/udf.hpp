#pragma once
#include "config.h"
#include "timarch.h"

inline void DEFINE_BOUNDARY();
inline void BEFORE_CONFIG();

void BEFORE_CONFIG(){
    cc::meshpath = "mesh/naca0012.txt";
    cc::testpath = "test.txt";
    cc::fieldpath = "field";
    cc::max_step = 200000;
    fatime::CFL = 1.0;
    fatime::USE_GLOBAL_DT = false;   // 定常: 当地时间步长
}

void DEFINE_BOUNDARY(){
    cc::VIL_DEFINE = {103.77, 9.08, 300.0, 3.3e6};
    cc::FAR_DEFINE = {103.77, 9.08, 300.0, 3.3e6};
    cc::POL_DEFINE = {3.3e6, 300.0};
}
