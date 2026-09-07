#pragma once
#include "mesh.h"

// 连接网格，计算单元面积、中心、面方向及固定几何缓存。
void initialize_geometry();
// 到最近壁面中点的距离 d；保存 1/d²，供 SA 源项和时间步使用。
void cache_wall_distance(cfd::Cell& cell);
