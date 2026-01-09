#include <iostream>
#include "src/ecmc/ecmc_frozen.h"
#include "src/mpi/MpiTopology.h"
#include <chrono>

#include "src/gauge/LocalGaugeField.h"
#include "src/mpi/HalosExchange.h"


void in_main_ecmc_frozen(char* argv[]) {
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
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc,&argv);
    mpi::MpiTopology topo(2);
    std::cout << "Rank :" << topo.rank << ", local :" << topo.local_rank << std::endl;
    mpi::GeometryFrozenMPI geo(4);
    mpi::LocalGaugeField field(4);
    mpi::shift::n_full_shifts(field, geo, 1, 1, topo);
    MPI_Finalize();
    return 0;
}