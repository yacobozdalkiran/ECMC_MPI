//
// Created by ozdalkiran-l on 1/21/26.
//

#include <iostream>
#include <chrono>

#include "../heatbath/heatbath.h"
#include "../geometry/Geometry.h"
#include "../observables/observables.h"
#include "../io/io.h"

//Prints the parameters of the run
void print_params(const RunParamsHb &rp) {
    std::cout << "==========================================" << std::endl;
    std::cout << "Lattice size : " << rp.L << "x" << rp.L << "x" << rp.L << "x" << rp.T << "\n";
    std::cout << "Start : " << (rp.cold_start ? "Cold" : "Hot") << "\n";
    std::cout << "Beta : " << rp.hp.beta << "\n";
    std::cout << "Number of samples : " << rp.hp.N_samples << "\n";
    std::cout << "Number of sweeps per sample : " << rp.hp.N_sweeps<< "\n";
    std::cout << "Number of hits per sweep : " << rp.hp.N_hits << "\n";
    std::cout << "Seed : " << rp.seed << "\n";
    std::cout << "==========================================" << std::endl;
}

//Write the plaquette in a txt file
void write_output(const std::vector<double> &meas, const RunParamsHb &rp) {
    int precision_filename = 1;
    std::string filename = "H_L"+std::to_string(rp.L)
                            + "T"+std::to_string(rp.T)
                            + "b"+io::format_double(rp.hp.beta, precision_filename)
                            + "c" + std::to_string(rp.cold_start)
                            + "Nsm" + std::to_string(rp.hp.N_sweeps)
                            + "Nh" + std::to_string(rp.hp.N_hits)
                            + "Ns" + std::to_string(rp.hp.N_samples);
    int precision = 10;
    io::save_double(meas, filename, precision);
}

//Reads the parameters of the input file into params
void read_params(RunParamsHb &params, const std::string &input) {
    try {
        io::load_params_hb(input, params);
    } catch (const std::exception& e) {
        std::cerr << "Error reading input : " << e.what() << std::endl;
    }
}

//Generates the samples using Heatbath
void generate_heatbath(const RunParamsHb &rp) {
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
    //Heatbath
    std::vector<double> meas = heatbath::samples(field, geo, rp.hp, rng);

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
    RunParamsHb params;
    read_params(params, argv[1]);
    print_params(params);

    auto start = std::chrono::high_resolution_clock::now();
    generate_heatbath(params);
    auto end = std::chrono::high_resolution_clock::now();
    long elapsed = std::chrono::duration_cast<std::chrono::seconds>(end-start).count();
    print_time(elapsed);

    return 0;
}