//
// Created by ozdalkiran-l on 1/16/26.
//

#include <iostream>
#include <chrono>

#include "../ecmc/ecmc.h"
#include "../gauge/GaugeField.h"
#include "../io/io.h"
#include "../observables/observables.h"

//Prints the params of rp
void print_params(const RunParamsSC &rp) {
    std::cout << "==========================================" << std::endl;
    std::cout << "Lattice size : " << rp.L << "x" << rp.L << "x" << rp.L << "x" << rp.T << "\n";
    std::cout << "Start : " << (rp.cold_start ? "Cold" : "Hot") << "\n";
    std::cout << "Beta : " << rp.ecmc_params.beta << "\n";
    std::cout << "Number of samples : " << rp.ecmc_params.N_samples << "\n";
    std::cout << "Theta sample : " << rp.ecmc_params.param_theta_sample << "\n";
    std::cout << "Theta refresh : " << rp.ecmc_params.param_theta_refresh << "\n";
    std::cout << "Loi de poisson : " << (rp.ecmc_params.poisson ? "Yes" : "No") << "\n";
    std::cout << "Epsilon set : " << rp.ecmc_params.epsilon_set << "\n";
    std::cout << "==========================================" << std::endl;
}

//Write the output
void write_output(const std::vector<double> &meas, const RunParamsSC &rp) {
    int precision_filename = 1;
    std::string filename = "L"+std::to_string(rp.L)
                            + "T"+std::to_string(rp.T)
                            + "b"+io::format_double(rp.ecmc_params.beta, precision_filename)
                            + "c" + std::to_string(rp.cold_start)
                            + "ts" + io::format_double(rp.ecmc_params.param_theta_sample, 0)
                            + "tr" + io::format_double(rp.ecmc_params.param_theta_refresh, 0)
                            + "Ns" + std::to_string(rp.ecmc_params.N_samples);
    int precision = 10;
    io::save_double(meas, filename, precision);
}

//Reads the parameters of the input file into params
void read_params(RunParamsSC &params, const std::string &input) {
    try {
        io::load_params_sc(input, params);
    } catch (const std::exception& e) {
        std::cerr << "Error reading input : " << e.what() << std::endl;
    }
}

//Prints the elapsed time
void print_time(long elapsed) {
    std::cout << "==========================================" << std::endl;
    std::cout << "Elapsed time : " << elapsed << "s\n";
    std::cout << "==========================================" << std::endl;
}

//Generates the samples using ECMC
void generate_sc(const RunParamsSC &rp) {
    int L = rp.L;
    int T = rp.T;

    //Random
    std::random_device rd;
    std::mt19937_64 rng(rd());

    //Initialization
    Geometry geo(L, T);
    GaugeField field(geo);
    if (!rp.cold_start) field.hot_start(rng);

    //ECMC
    std::vector<double> meas = ecmc::samples_improved(field, geo, rp.ecmc_params, rng);

    //Output
    write_output(meas, rp);
}

//Binary gauge_single_core to generate single core lattice gauge configurations with ECMC
int in_main_gauge_single_core(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file.txt>" << std::endl;
        return 1;
    }

    RunParamsSC params;
    read_params(params, argv[1]);
    print_params(params);
    auto start = std::chrono::high_resolution_clock::now();
    generate_sc(params);
    auto end = std::chrono::high_resolution_clock::now();
    long elapsed = std::chrono::duration_cast<std::chrono::seconds>(end-start).count();
    print_time(elapsed);
    return 0;
}

//Checks the gauge invariance of the plaquette
void check_plaquette(const RunParamsSC &rp) {
    int L = rp.L;
    int T = rp.T;

    //Random
    std::random_device rd;
    std::mt19937_64 rng(rd());

    //Initialization
    Geometry geo(L, T);
    GaugeField field(geo);
    if (!rp.cold_start) field.hot_start(rng);

    std::cout << "Plaquette = " << observables::mean_plaquette(field, geo) << '\n';
    std::cout << "Gauge transform...\n";
    observables::gauge_transform(field, geo, rng);
    std::cout << "Plaquette = " << observables::mean_plaquette(field, geo) << '\n';
}

//Binary to check plaquette with gauge invariance
int in_main_check_plaquette_gauge(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file.txt>" << std::endl;
        return 1;
    }

    RunParamsSC params;
    read_params(params, argv[1]);
    print_params(params);
    auto start = std::chrono::high_resolution_clock::now();
    check_plaquette(params);
    auto end = std::chrono::high_resolution_clock::now();
    long elapsed = std::chrono::duration_cast<std::chrono::seconds>(end-start).count();
    print_time(elapsed);
    return 0;
}

int main(int argc, char* argv[]) {
    int status = in_main_gauge_single_core(argc, argv);
    return status;
}