#include "config.h"
#include "timarch.h"

inline void DEFINE_BOUNDARY();

inline void BEFORE_CONFIG();

// 前期设置
void BEFORE_CONFIG(){
    cc::meshpath = "mesh/tunnel.txt";
    cc::testpath = "test.txt";
    cc::turbulence = "SA";
    cc::max_step = 100;
    fatime::CFL = 0.5;   // 时间步长安全系数(可在此调整, 更小更稳定)
    fatime::USE_GLOBAL_DT = true;   // 瞬态模式: 全局最小时间步推进
}

// 定义边界条件
void DEFINE_BOUNDARY(){
    cc::SA tur_POL = {};
    tur_POL.miubl = 1e-6;
    cc::SA tur_VIL = {};
    tur_VIL.miubl = 1e-6;
    cc::POL_DEFINE = {90000.0,300,tur_POL};   // 出口压力90000, 顺压梯度
    cc::VIL_DEFINE = {10,0,300,101325,tur_VIL};
}