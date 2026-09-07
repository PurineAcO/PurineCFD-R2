#pragma once

namespace parallel {
// Use visible physical cores by default; OMP_NUM_THREADS remains authoritative.
const char* configure_threads();
} // namespace parallel
