#pragma once
#include "classconfig.h"

// 找到全部邻接节点
void findnode(cc::cell_class& cell);
// 计算网格体积
double volume(cc::cell_class& cell);
// 计算网格质心
void center(cc::cell_class& cell);
// 几何分析主程序
void geometrymain();