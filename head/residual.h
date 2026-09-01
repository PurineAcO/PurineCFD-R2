#pragma once
#include <vector>

namespace res {

    // 累加 5 个方程的残差平方
    void accumulate(double resSq[5], const std::vector<double>& conser_RK);

    // 报告残差
    void report(double resSq[5], int step, double init[5]);

}
