
#include <chrono>
#include <iomanip>
#include <iostream>

#include "../ecmc/ecmc_mpi_cb.h"
#include "../gauge/GaugeField.h"
#include "../io/io.h"
#include "../mpi/Halo.h"
#include "../mpi/HalosCB.h"
#include "../mpi/HalosExchange.h"
#include "../mpi/MpiTopology.h"

// Generates samples of gauge confs with 4d shifts
void generate_fullshift(const RunParamsCB& run_params) {
    int L = run_params.L_core;
    int n_core_dims = run_params.n_core_dims;

    // Create the topology
    mpi::MpiTopology topo(n_core_dims);

    // Rng
    // std::random_device rd;
    std::mt19937_64 rng(run_params.seed + topo.rank);

    // Create the geometry and the field (cold)
    GeometryCB geo(L);
    GaugeField field(geo);

    // Initialize the field
    if (!run_params.cold_start) field.hot_start(rng);

    // Create the ECMC halos
    HalosCB halocb(geo);
    mpi::ecmccb::fill_and_exchange(field, geo, halocb, topo);
    // Create the observables halos
    HaloObs halo_obs(geo);

    // Shift params and halo creation
    int L_shift = L / 2;  // Max shift length
    Halo halo_shift(L_shift, geo);
    ShiftParams sp;  // pos and coord and L_shift unset by default
    if (run_params.stype_pos)
        sp.stype = pos;
    else
        sp.stype = neg;
    sp.L_shift = L_shift;
    // For random generation of shifts
    std::uniform_int_distribution<int> dist_shift(1, L_shift);
    int xshift{}, yshift{}, zshift{}, tshift{};

    // Global shift parameters
    int n_shift = run_params.n_shift;

    // ECMC params
    ECMCParams ecmc_params = run_params.ecmc_params;
    int N_samples = ecmc_params.N_samples;

    // Space reserved for results
    std::vector<std::vector<double>> plaquette(2 * n_shift, std::vector<double>(N_samples, 0.0));

    // Print parameters
    if (topo.rank == 0) {
        std::cout << "==========================================" << std::endl;
        std::cout << "Total lattice size : " << L * n_core_dims << "^4\n";
        std::cout << "Local lattice size : " << L << "^4\n";
        std::cout << "Beta : " << ecmc_params.beta << "\n";
        std::cout << "Total number of shifts : " << n_shift << "\n";
        std::cout << "Number of samples per shift : " << ecmc_params.N_samples << "\n";
        std::cout << "Total number of samples : " << 2 * ecmc_params.N_samples * n_shift << "\n";
        std::cout << "Seed : " << run_params.seed << "\n";
        std::cout << "==========================================" << std::endl;
    }

    for (int gshiftc = 0; gshiftc < n_shift; gshiftc++) {
        // Random generation of shifts
        if (topo.rank == 0) {
            xshift = dist_shift(rng);
            yshift = dist_shift(rng);
            zshift = dist_shift(rng);
            tshift = dist_shift(rng);
        }

        // Synchronization
        MPI_Bcast(&xshift, sizeof(int), MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&yshift, sizeof(int), MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&zshift, sizeof(int), MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&tshift, sizeof(int), MPI_INT, 0, MPI_COMM_WORLD);

        // Shifts
        if (topo.rank == 0) std::cout << "Shift " << gshiftc << " (X)";
        sp.coord = X;
        sp.L_shift = xshift;
        mpi::shiftcb::shift(field, geo, halo_shift, topo, sp);
        if (topo.rank == 0) std::cout << "(Y)";
        sp.coord = Y;
        sp.L_shift = yshift;
        mpi::shiftcb::shift(field, geo, halo_shift, topo, sp);
        if (topo.rank == 0) std::cout << "(Z)";
        sp.coord = Z;
        sp.L_shift = zshift;
        mpi::shiftcb::shift(field, geo, halo_shift, topo, sp);
        if (topo.rank == 0) std::cout << "(T)\n";
        sp.coord = T;
        sp.L_shift = tshift;
        mpi::shiftcb::shift(field, geo, halo_shift, topo, sp);

        // Fill the halo for ECMC
        mpi::ecmccb::fill_and_exchange(field, geo, halocb, topo);
        parity active_parity = even;
        plaquette[2 * gshiftc] = mpi::ecmccb::samples_improved(field, geo, ecmc_params, rng,
                                                               halo_obs, topo, active_parity);
        active_parity = odd;
        mpi::ecmccb::fill_and_exchange(field, geo, halocb, topo);
        plaquette[2 * gshiftc + 1] = mpi::ecmccb::samples_improved(field, geo, ecmc_params, rng,
                                                                   halo_obs, topo, active_parity);

        if (topo.rank == 0) {
            // Flatten the plaquette vector
            std::vector<double> plaquette_flat(2 * n_shift * N_samples);
            for (int i = 0; i < 2 * n_shift; i++) {
                for (int j = 0; j < N_samples; j++) {
                    plaquette_flat[i * N_samples + j] = plaquette[i][j];
                }
            }
            // Write the output
            int precision_filename = 1;
            std::string filename =
                "EMCB_" + std::to_string(L * n_core_dims) + "b" +
                io::format_double(ecmc_params.beta, precision_filename) + "Ns" +
                std::to_string(run_params.n_shift) + "c" + std::to_string(run_params.cold_start) +
                "ts" +
                io::format_double(run_params.ecmc_params.param_theta_sample, precision_filename) +
                "tr" +
                io::format_double(run_params.ecmc_params.param_theta_refresh, precision_filename);
            int precision = 10;
            io::save_double(plaquette_flat, filename, precision);
        }
    }
}

// Reads the parameters of input file into RunParams struct
void read_params(RunParamsCB& params, int rank, const std::string& input) {
    if (rank == 0) {
        try {
            io::load_params(input, params);
        } catch (const std::exception& e) {
            std::cerr << "Error reading input : " << e.what() << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    // Synchronizing input parameters accross all nodes
    MPI_Bcast(&params, sizeof(RunParamsCB), MPI_BYTE, 0, MPI_COMM_WORLD);
}

// Main function
int main(int argc, char* argv[]) {
    // Initialization and check args
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // Trying to read input
    if (argc < 2) {
        if (rank == 0) {
            std::cerr << "Usage: " << argv[0] << " <input_file.txt>" << std::endl;
        }
        MPI_Finalize();
        return 1;
    }

    // Charging the parameters of the run
    RunParamsCB params;
    read_params(params, rank, argv[1]);

    // Measuring time
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
