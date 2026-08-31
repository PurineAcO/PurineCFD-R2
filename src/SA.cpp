// #include "classconfig.h"
#include "SA.h"

// 计算湍流粘度比chi
double SA::chi(double miubl, double miu,double rho){return (rho*miubl)/miu; }

// 计算粘度阻尼函数fv1
double SA::fv1(double chi){return (chi*chi*chi)/(chi*chi*chi+SA::Cv1*SA::Cv1*SA::Cv1);}

// S-A主函数
// double SA::SA_main(cc::cell_class &cell){

// }