#include "residual.h"
#include <cmath>
#include <cstdio>

namespace res {

void accumulate(double resSq[5], const std::vector<double>& conser_RK){
    for(int s = 0; s < 5; s++){
        resSq[s] += conser_RK[s] * conser_RK[s];
    }
}

void report(double resSq[5], int step, double init[5]){
    static bool header = false;
    if(!header){
        printf("%6s %12s %12s %12s %12s %12s\n",
               "step", "continuous", "vx", "vy", "energy", "nut");
        header = true;
    }
    double L2[5];
    for(int s = 0; s < 5; s++){ L2[s] = std::sqrt(resSq[s]); }
    if(step == 1){ for(int s = 0; s < 5; s++){ init[s] = L2[s]; } }
    printf("%6d", step);
    for(int s = 0; s < 5; s++){
        double norm = (init[s] > 0) ? L2[s] / init[s] : 0.0;
        printf(" %12.6e", norm);
    }
    printf("\n");
}

}
