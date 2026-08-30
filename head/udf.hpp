#include "config.h"

inline void DEFINE_BOUNDARY();

// 定义边界条件
void DEFINE_BOUNDARY(){
    cc::SA tur_POL = {};
    tur_POL.miubl = 1e-4;
    cc::SA tur_VIL = {};
    tur_VIL.miubl = 1e-4;
    cc::POL_DEFINE = {101325.0,300,tur_POL};
    cc::VIL_DEFINE = {10,0,300,101325,tur_VIL};
}