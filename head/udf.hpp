#pragma once
#include "config.h"
#include "timarch.h"

inline void DEFINE_BOUNDARY();
inline void BEFORE_CONFIG();

void BEFORE_CONFIG(){
    cc::meshpath = "mesh/tunnel_uni.txt";
    cc::testpath = "test.txt";
    cc::fieldpath = "field";
    cc::max_step = 1000000;
    fatime::CFL = 2.0;
    fatime::USE_GLOBAL_DT = true;
}

void DEFINE_BOUNDARY(){
    cc::POL_DEFINE = {90000.0, 300.0};
    cc::VIL_DEFINE = {100.0, 0.0, 300.0, 101325.0};
}
