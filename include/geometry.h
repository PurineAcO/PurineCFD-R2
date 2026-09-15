#pragma once
#include "classconfig.h"

// 找到全部邻接节点
void findnode(cc::cell_class& cell);
// 计算网格体积(四个三角形面积和的一半)
void volume(cc::cell_class& cell);
// 计算网格质心(按三角形面积加权)
void center(cc::cell_class& cell);
// 计算壁面距离
void sad(cc::cell_class& cell);
// 几何分析主程序
bool geometrymain();
