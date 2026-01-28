#include <iostream>
#include <chrono>
#include <iomanip>

#include "../heatbath/heatbath_mpi.h"
#include "../gauge/GaugeField.h"
#include "../io/io.h"
#include "../mpi/Halo.h"
#include "../mpi/HalosExchange.h"
#include "../mpi/MpiTopology.h"


//Generates samples of gauge confs with 4d shifts
void generate_fullshift(const RunParamsHbMPI &run_params) {
    int L = run_params.L_core;
    int n_core_dims = run_params.n_core_dims;

    //Create the topology
    mpi::MpiTopology topo(n_core_dims);

    //Rng
    //std::random_device rd;
    std::mt19937_64 rng(run_params.seed+topo.rank);

    //Create the geometry and the field (cold)
    GeometryHaloECMC geo(L);
    GaugeField field(geo);

    //Initialize the field
    if (!run_params.cold_start) field.hot_start(rng);

    //Create the ECMC halos
    HaloECMC halo_ecmc(geo);
    mpi::ecmc::fill_and_exchange(field, geo, halo_ecmc, topo);
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

    //Hb params
    HbParams hb_params = run_params.hp;
    int N_samples = hb_params.N_samples;

    //Space reserved for results
    std::vector<std::vector<double>> plaquette(n_shift);

    //Print parameters
    if (topo.rank == 0) {
        std::cout << "==========================================" << std::endl;
        std::cout << "Total lattice size : " << L*n_core_dims  << "^4\n";
        std::cout << "Local lattice size : " << L << "^4\n";
        std::cout << "Beta : " << hb_params.beta << "\n";
        std::cout << "Total number of shifts : " << n_shift << "\n";
        std::cout << "Offset of shifts : " << run_params.L_shift << "\n";
        std::cout << "Number of samples per shift : " << hb_params.N_samples << "\n";
        std::cout << "Total number of samples : " << hb_params.N_samples*n_shift << "\n";
        std::cout << "Seed : " << run_params.seed << "\n";
        std::cout << "==========================================" << std::endl;
        int g = std::gcd(run_params.L_core, run_params.L_shift);
        if (not((L>4*g))) std::cout << "/!\\ WARNING : Not ergodic ! Change L_core or L_shift /!\\ \n";
    }

    for (int gshiftc = 0; gshiftc<n_shift; gshiftc++) {
        //Shifts
        if (topo.rank == 0) std::cout << "Shift " << gshiftc << " (X)";
        sp.coord=X;
        mpi::shift::shift(field, geo, halo_shift, topo, sp);
        if (topo.rank == 0) std::cout << "(Y)";
        sp.coord=Y;
        mpi::shift::shift(field, geo, halo_shift, topo, sp);
        if (topo.rank == 0) std::cout << "(Z)";
        sp.coord=Z;
        mpi::shift::shift(field, geo, halo_shift, topo, sp);
        if (topo.rank == 0) std::cout <<  "(T)\n";
        sp.coord=T;
        mpi::shift::shift(field, geo, halo_shift, topo, sp);

        //Fill the halo for ECMC
        mpi::ecmc::fill_and_exchange(field, geo, halo_ecmc, topo);

        plaquette[gshiftc] = mpi::heatbath::samples(field, geo, halo_obs, topo, hb_params, rng);
    }

    if (topo.rank == 0) {
        //Flatten the plaquette vector
        std::vector<double> plaquette_flat(n_shift*N_samples);
        for (int i = 0; i<n_shift; i++) {
            for (int j = 0; j<N_samples; j++) {
                plaquette_flat[i*N_samples + j] = plaquette[i][j];
            }
        }
        //Write the output
        int precision_filename = 1;
        std::string filename = "HBM_L" + std::to_string(L*n_core_dims)
                                + "b" + io::format_double(hb_params.beta, precision_filename)
                                + "Ls" + std::to_string(run_params.L_shift)
                                + "Ns" + std::to_string(run_params.n_shift)
                                + "c" + std::to_string(run_params.cold_start)
                                + "Nsm" + std::to_string(hb_params.N_sweeps)
                                + "Nh" + std::to_string(hb_params.N_hits)
                                + "Ns" + std::to_string(hb_params.N_samples*run_params.n_shift);
        int precision = 10;
        io::save_double(plaquette_flat, filename, precision);
    }
}

//Reads the parameters of input file into RunParams struct
void read_params(RunParamsHbMPI &params, int rank, const std::string &input) {
    if (rank == 0) {
        try {
            io::load_params(input, params);
        } catch (const std::exception& e) {
            std::cerr << "Error reading input : " << e.what() << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    //Synchronizing input parameters accross all nodes
    MPI_Bcast(&params, sizeof(RunParamsHbMPI), MPI_BYTE, 0, MPI_COMM_WORLD);
}

//Main function
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
    RunParamsHbMPI params;
    read_params(params, rank, argv[1]);
    //Measuring time
    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    generate_fullshift(params);

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
