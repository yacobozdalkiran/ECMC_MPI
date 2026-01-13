#include <iostream>
#include <chrono>
#include "src/ecmc/ecmc_frozen.h"
#include "src/ecmc/ecmc_mpi.h"
#include "src/gauge/GaugeField.h"
#include "src/mpi/HalosExchange.h"
#include "src/mpi/MpiTopology.h"


void in_main_ecmc_frozen(char* argv[]) {
    int L = atoi(argv[1]);
    int T = atoi(argv[1]);
    mpi::GeometryFrozen geo(L);
    GaugeField field(geo);
    std::random_device rd;
    std::mt19937_64 rng(rd());
    auto start = std::chrono::high_resolution_clock::now();
    ecmc_frozen::samples_improved(field, geo, atoi(argv[2]), atoi(argv[3]), atoi(argv[4]),atoi(argv[5]), atoi(argv[6]),atoi(argv[7]),rng);
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = end - start;
    std::cout << "Elapsed time : " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count() << " ms" << std::endl;
}

void in_main_MPI_test(int argc, char* argv[]) {
    MPI_Init(&argc,&argv);


    int L = 6;
    mpi::MpiTopology topo(2);
    mpi::GeometryFrozen geo(L);
    GaugeField field(geo);
    //field.hot_start(rng);
    HaloObs halo_obs(geo);

    //ECMC params
    double beta = 6.0;
    int N_samples = 50;
    double param_theta_sample = 1000;
    double param_theta_refresh = 300;
    bool poisson = false;
    double epsilon_set = 0.15;

    std::mt19937_64 rng(123+topo.rank);

    //Shift params
    int L_shift = 2;
    halo_coord coord = X;
    shift_type stype = pos;
    Halo halo_shift(L_shift, geo, coord);

    mpi::ecmc::samples_improved(field, geo, beta, N_samples, param_theta_sample, param_theta_refresh, poisson, epsilon_set, rng, halo_obs, topo);
    for (int i =0; i<10; i++) {

        if (topo.rank == 0) std::cout << "Shifting...\n";
        coord = X;
        mpi::shift::shift(field, geo, halo_shift, topo, L_shift, coord, stype);
        mpi::ecmc::samples_improved(field, geo, beta, N_samples, param_theta_sample, param_theta_refresh, poisson, epsilon_set, rng, halo_obs, topo);

        if (topo.rank == 0) std::cout << "Shifting...\n";
        coord = Y;
        mpi::shift::shift(field, geo, halo_shift, topo, L_shift, coord, stype);
        mpi::ecmc::samples_improved(field, geo, beta, N_samples, param_theta_sample, param_theta_refresh, poisson, epsilon_set, rng, halo_obs, topo);

        if (topo.rank == 0) std::cout << "Shifting...\n";
        coord = Z;
        mpi::shift::shift(field, geo, halo_shift, topo, L_shift, coord, stype);
        mpi::ecmc::samples_improved(field, geo, beta, N_samples, param_theta_sample, param_theta_refresh, poisson, epsilon_set, rng, halo_obs, topo);

        if (topo.rank == 0) std::cout << "Shifting...\n";
        coord = T;
        mpi::shift::shift(field, geo, halo_shift, topo, L_shift, coord, stype);

        mpi::ecmc::samples_improved(field, geo, beta, N_samples, param_theta_sample, param_theta_refresh, poisson, epsilon_set, rng, halo_obs, topo);
    }
    MPI_Finalize();
}

int main(int argc, char* argv[]) {
    in_main_MPI_test(argc, argv);
    return 0;
}