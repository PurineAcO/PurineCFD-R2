#pragma once

// 建立日志目录并重定向stdout
bool open_log(const char* path);
// 输出流场
bool dump_field(int step);
