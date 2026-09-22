#pragma once



// 边界条件参数: 由config.json读入, 供boundary.cpp使用

// 远场自由流(亚声速)
struct far_condition{
    double u = 0.0;      // x方向速度
    double v = 0.0;      // y方向速度
    double T = 300.0;    // 温度
    double p = 101325.0; // 压力
};

inline far_condition FAR_DEFINE;
