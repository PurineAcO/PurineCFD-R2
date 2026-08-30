#pragma once

namespace cc {

    struct physics{

    // 重要变量
    double rho; // 密度
    double u;   // u速度
    double v;   // v速度
    double T;   // 温度

    // 衍生物
    double a;   // 声速
    double mu;  // 空气粘度
    double p;   // 压强

    };// 物理量

    struct SA{
        double miubl; // 工作变量miubl
    };// S-A湍流变量

    using TUR = SA ;  // 在此修改湍流模型

    // 求解器设置
    inline int cell_num = 0;                          // 网格总数
    inline int face_num = 0;                          // 面总数
    inline int type_total = 0;                        // 边界条件类型总数
    inline const int HALO = 2;                        // HALO网格层数(可能弃用)
    inline const char* meshpath = "mesh/tunnel.txt";  // 网格文件位置
    inline const char* testpath = "test.txt";         // 测试文件位置

    // 网格和边界条件约定
    inline const short INTER = 0;                     // 内部面
    inline const short WALL = 1;                      // 无滑移壁面
    inline const short POL = 2;                       // 压力出口
    inline const short VIL = 3;                       // 速度入口
    inline const short FAR = 4;                       // 压力远场

    // 临时的几何参数
    inline const double H = 2.4;                      // 方腔高度

    // 流动物理条件
    inline struct VIL_condition{double u;double v;            // 速度分量
                                double T;                     // 入口温度,必须是静温
                                double p;                     // 来流静压,亚声速不适用
                                TUR tur;                      // 湍流模型
                                }VIL_DEFINE;   // 速度入口
    inline struct POL_condition{double p;                    // 出口静压
                                double T;                    // 出口温度,必须是总温
                                TUR tur;                      // 湍流模型
                                }POL_DEFINE;   // 压力出口

    // 物理学常数
    inline const double gamma = 1.4;
    inline const double R = 287.05;
    inline const double mu_ref = 1.716e-5;
    inline const double T_ref = 273.15;
    inline const double T_s = 110.4;


}