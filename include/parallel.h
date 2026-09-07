#pragma once

namespace parallel {
// 优先采用 OMP_NUM_THREADS；未设置时使用进程可用的物理核心数。
const char* configure_threads();
} // namespace parallel
