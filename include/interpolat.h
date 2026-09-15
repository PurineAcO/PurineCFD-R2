#pragma once
#include "classconfig.h"

// 面上中心差分插值; 内部面只由nei[0]侧单元调用
void interpolate_mid(cc::cell_class& cell);
