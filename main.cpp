#include <iostream>
#include "src/ecmc/ecmc_frozen.h"
#include <chrono>

int main(int argc, char* argv[]) {
    int L = atoi(argv[1]);
    int T = atoi(argv[1]);
    GaugeField field(L, T);
    GeometryFrozen geo(L, T);
    std::random_device rd;
    std::mt19937_64 rng(rd());
    auto start = std::chrono::high_resolution_clock::now();
    ecmc_frozen::samples_improved(field, geo, atoi(argv[2]), atoi(argv[3]), atoi(argv[4]),atoi(argv[5]), atoi(argv[6]),atoi(argv[7]),rng);
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = end - start;
    std::cout << "Elapsed time : " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count() << " ms" << std::endl;
    return 0;
}