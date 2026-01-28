//
// Created by ozdalkiran-l on 1/20/26.
//
#include "../metropolis/metropolis.h"
#include "../gauge/GaugeField.h"
#include "../observables/observables.h"
#include <chrono>
#include "../io/io.h"
#include <iostream>

//Prints the params of rp
void print_params(const RunParamsMetro &rp) {
    std::cout << "==========================================" << std::endl;
    std::cout << "Lattice size : " << rp.L << "x" << rp.L << "x" << rp.L << "x" << rp.T << "\n";
    std::cout << "Start : " << (rp.cold_start ? "Cold" : "Hot") << "\n";
    std::cout << "Beta : " << rp.mp.beta << "\n";
    std::cout << "Number of samples : " << rp.mp.N_samples << "\n";
    std::cout << "Number of sweeps per sample : " << rp.mp.N_sweeps_meas << "\n";
    std::cout << "Number of hits per sweep : " << rp.mp.N_hits << "\n";
    std::cout << "Epsilon : " << rp.mp.epsilon << "\n";
    std::cout << "Seed : " << rp.seed << "\n";
    std::cout << "==========================================" << std::endl;
}

//Write the output
void write_output(const std::vector<double> &meas, const RunParamsMetro &rp) {
    int precision_filename = 1;
    std::string filename = "M_L"+std::to_string(rp.L)
                            + "T"+std::to_string(rp.T)
                            + "b"+io::format_double(rp.mp.beta, precision_filename)
                            + "c" + std::to_string(rp.cold_start)
                            + "Nh" + std::to_string(rp.mp.N_hits)
                            + "Nsm" + std::to_string(rp.mp.N_sweeps_meas)
                            + "Ns" + std::to_string(rp.mp.N_samples);
    int precision = 10;
    io::save_double(meas, filename, precision);
}

//Reads the parameters of the input file into params
void read_params(RunParamsMetro &params, const std::string &input) {
    try {
        io::load_params(input, params);
    } catch (const std::exception& e) {
        std::cerr << "Error reading input : " << e.what() << std::endl;
    }
}

//Generates the samples using Metropolis
void generate_metro(const RunParamsMetro &rp) {
    int L = rp.L;
    int T = rp.T;

    //Random
    //std::random_device rd;
    std::mt19937_64 rng(rp.seed);

    //Initialization
    Geometry geo(L, T);
    GaugeField field(geo);
    if (!rp.cold_start) field.hot_start(rng);

    std::cout << "Max drift from unitarity : " << observables::max_drift_det(field, geo) << "\n";
    //Metropolis
    std::vector<double> meas = metropolis::samples(field, geo, rp.mp, rng);

    //Check SU3
    std::cout << "Max drift from unitarity : " << observables::max_drift_det(field, geo) << "\n";

    //Output
    write_output(meas, rp);
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file.txt>" << std::endl;
        return 1;
    }
    RunParamsMetro params;
    read_params(params, argv[1]);
    print_params(params);
    auto start = std::chrono::high_resolution_clock::now();
    generate_metro(params);
    auto end = std::chrono::high_resolution_clock::now();
    long elapsed = std::chrono::duration_cast<std::chrono::seconds>(end-start).count();
    print_time(elapsed);
    return 0;
}