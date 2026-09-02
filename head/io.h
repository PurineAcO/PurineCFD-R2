#pragma once

// 输出全场控制体
void dump_field(double t);

// 保存checkpoint
void save_checkpoint(double t);

// 读取checkpoint;成功返回true并恢复全场状态
bool load_checkpoint(double& t);
