#pragma once

namespace res {
    void report_update(int step);
    double current_residual();       // 原残差: 密度场更新RMS
    double current_residual_all();   // 新残差: 各胞四量差平方和均值
}
