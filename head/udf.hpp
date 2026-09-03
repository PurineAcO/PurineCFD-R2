#pragma once
#include "config.h"
#include "timarch.h"

inline void DEFINE_BOUNDARY();
inline void BEFORE_CONFIG();

void BEFORE_CONFIG(){
    cc::meshpath = "mesh/cyl.txt";
    cc::testpath = "test.txt";
    cc::fieldpath = "field";
    cc::max_step = 200000;
    fatime::CFL = 1.0;
    fatime::USE_GLOBAL_DT = false;   // 定常: 当地时间步长
}

void DEFINE_BOUNDARY(){
    cc::VIL_DEFINE = {520.83, 0.0, 300.0, 101325.0};  // Ma=1.5 超声速自由流
    cc::FAR_DEFINE = {520.83, 0.0, 300.0, 101325.0};
    cc::POL_DEFINE = {101325.0, 300.0};
}
