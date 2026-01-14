//
// Created by ozdalkiran-l on 1/14/26.
//

#include <iostream>
#include <chrono>

#include "../ecmc/ecmc_mpi.h"
#include "../gauge/GaugeField.h"
#include "../io/io.h"
#include "../mpi/Halo.h"
#include "../mpi/HalosExchange.h"
#include "../mpi/MpiTopology.h"

void generate(const RunParams &run_params) {
    int L = run_params.L_core;
    int n_core_dims = run_params.n_core_dims;

    //Create the topology
    mpi::MpiTopology topo(n_core_dims);

    //Rng
    std::random_device rd;
    std::mt19937_64 rng(rd());

    //Create the geometry and the field (cold)
    mpi::GeometryFrozen geo(L);
    GaugeField field(geo);

    //Initialize the field
    if (!run_params.cold_start) field.hot_start(rng);

    //Create the observables halos
    HaloObs halo_obs(geo);

    //Shift params and halo creation
    int L_shift = run_params.L_shift; //Max shift length
    Halo halo_shift(L_shift, geo);
    ShiftParams sp; //pos and coord and L_shift unset by default
    if (run_params.stype_pos) sp.stype = pos;
    else sp.stype = neg;
    sp.L_shift = L_shift;

    //Global shift parameters
    int n_shift = run_params.n_shift;
    int total_shifts = n_shift*4;

    //ECMC params
    ECMCParams ecmc_params = run_params.ecmc_params;
    int N_samples = ecmc_params.N_samples;

    //Space reserved for results
    std::vector<std::vector<double>> plaquette(total_shifts);

    //Print parameters
    if (topo.rank == 0) {
        std::cout << "==========================================" << std::endl;
        std::cout << "Total lattice size : " << L*n_core_dims  << "\n";
        std::cout << "Local lattice size : " << L << "\n";
        std::cout << "Beta : " << ecmc_params.beta << "\n";
        std::cout << "Total number of shifts : " << n_shift*4 << "\n";
        std::cout << "Number of samples per shift : " << ecmc_params.N_samples << "\n";
        std::cout << "==========================================" << std::endl;
    }

    for (int gshiftc = 0; gshiftc<n_shift; gshiftc++) {
        if (topo.rank == 0) std::cout << "Shift " << gshiftc*4 << " (X)...\n";
        sp.coord=X;
        mpi::shift::shift(field, geo, halo_shift, topo, sp);
        plaquette[gshiftc*4] = mpi::ecmc::samples_improved(field, geo, ecmc_params, rng, halo_obs, topo);

        if (topo.rank == 0) std::cout << "Shift " << gshiftc*4+1 << " (Y)...\n";
        sp.coord=Y;
        mpi::shift::shift(field, geo, halo_shift, topo, sp);
        plaquette[gshiftc*4+1] = mpi::ecmc::samples_improved(field, geo, ecmc_params, rng, halo_obs, topo);

        if (topo.rank == 0) std::cout << "Shift " << gshiftc*4+2 << " (Z)...\n";
        sp.coord=Z;
        mpi::shift::shift(field, geo, halo_shift, topo, sp);
        plaquette[gshiftc*4+2] = mpi::ecmc::samples_improved(field, geo, ecmc_params, rng, halo_obs, topo);

        if (topo.rank == 0) std::cout << "Shift " << gshiftc*4+3 << " (T)...\n";
        sp.coord=T;
        mpi::shift::shift(field, geo, halo_shift, topo, sp);
        plaquette[gshiftc*4+3] = mpi::ecmc::samples_improved(field, geo, ecmc_params, rng, halo_obs, topo);
    }

    if (topo.rank == 0) {
        //Flatten the plaquette vector
        std::vector<double> plaquette_flat(total_shifts*N_samples);
        for (int i = 0; i<total_shifts; i++) {
            for (int j = 0; j<N_samples; j++) {
                plaquette_flat[i*N_samples + j] = plaquette[i][j];
            }
        }
        //Write the output
        std::string filename = "L"+std::to_string(L*n_core_dims)
                                + "b"+std::to_string(ecmc_params.beta)
                                + "Ls" + std::to_string(run_params.L_shift)
                                + "Ns" + std::to_string(run_params.n_shift)
                                + "dir" + std::to_string(run_params.stype_pos);
        int precision = 10;
        io::save_double_params(plaquette_flat, run_params, filename, precision);
    }
}

void read_params(RunParams &params, int rank, const std::string &input) {
    if (rank == 0) {
        try {
            io::load_params(input, params);
        } catch (const std::exception& e) {
            std::cerr << "Error reading input : " << e.what() << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    //Synchronizing input parameters accross all nodes
    MPI_Bcast(&params, sizeof(RunParams), MPI_BYTE, 0, MPI_COMM_WORLD);
}

int main(int argc, char* argv[]) {
    //Initialization and check args
    MPI_Init(&argc,&argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //Trying to read input
    if (argc < 2) {
        if (rank == 0) {
            std::cerr << "Usage: " << argv[0] << " <input_file.txt>" << std::endl;
        }
        MPI_Finalize();
        return 1;
    }

    //Charging the parameters of the run
    RunParams params;
    read_params(params, rank, argv[1]);

    //Measuring time
    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    generate(params);

    MPI_Barrier(MPI_COMM_WORLD);
    double end_time = MPI_Wtime();

    if (rank == 0) {
        double total_time = end_time - start_time;
        std::cout << std::fixed << std::setprecision(4);
        std::cout << "\n==========================================" << std::endl;
        std::cout << " Total execution time : " << total_time << " seconds" << std::endl;
        std::cout << "==========================================\n" << std::endl;
    }

    MPI_Finalize();
    return 0;
}
