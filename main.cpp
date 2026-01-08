#include <iostream>
#include "src/ecmc/EcmcUpdater.h"
#include <chrono>

int main() {
    int L = 4;
    int T = 4;
    GaugeField field(L, T);
    Geometry geo(L, T);
    EcmcUpdater ecmc(123);
    auto start = std::chrono::high_resolution_clock::now();
    ecmc.ecmc_samples_improved(field, geo, 3.0, 10000, 100,30,0,0.15);
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = end - start;
    std::cout << "Elapsed time : " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count() << " ms" << std::endl;
    return 0;
}