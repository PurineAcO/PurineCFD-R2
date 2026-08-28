#include "classconfig.h"

// 计算网格体积,使用四三角形法
double volume(cc::cell_class& cell);

// 计算网格质心,使用面积加权法
void center(cc::cell_class& cell);

// 计算到壁面的距离
void walldistance(cc::cell_class& cell);

// 计算到壁面的距离(方腔流动代替品)
void walldistance_alternative(cc::cell_class &cell);

void geometry();