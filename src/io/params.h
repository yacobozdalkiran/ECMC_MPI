//
// Created by ozdalkiran-l on 1/14/26.
//

#ifndef ECMC_MPI_PARAMS_H
#define ECMC_MPI_PARAMS_H

#include <string>

struct ECMCParams {
    double beta = 6.0;
    double param_theta_sample = 100;
    double param_theta_refresh = 50;
    bool poisson = false;
    double epsilon_set = 0.15;
};

struct HbParams {
    double beta = 6.0;
    int N_hits = 1;
    int N_sweeps = 1;
};

struct RunParamsECB {
    int L_core = 6;
    int n_core_dims = 2;
    bool cold_start = true;
    int seed = 123;
    int N_shift = 2;            // Number of shifts
    int N_ov_sweep = 2;               // Nb ov sweeps
    int N_ov_hit = 1;           // Nb ov hits
    ECMCParams ecmc_params{};   // Params of the ECMC for each even/odd update
    int N_shift_plaquette = 2;  // Measure plaquette every N_shift_plaquette_shift
    bool topo = true;
    int N_shift_topo = 10;  // Measure topo charge every N_shift_topo shift
    int N_steps_gf = 10;
    int N_rk_steps = 40;
    std::string run_name = "c";
    std::string run_dir = "data";
    int save_each_shifts = 2;  // save confs/measures/seed each
};

struct RunParamsHbCB {
    int L_core = 6;
    int n_core_dims = 2;
    bool cold_start = true;
    int N_shift = 2;  // Number of shifts
    int N_ov_sweep = 2;               // Nb ov sweeps
    int N_ov_hit = 1;           // Nb ov hits
    int seed = 123;
    HbParams hp{};  // Hb params for each even/odd update
    bool topo = true;
    int N_shift_plaquette = 2;  // Measure plaquette every N_shift_plaquette_shift
    int N_shift_topo = 10;      // Measure topo charge every N_shift_topo shift
    int N_steps_gf = 10;
    int N_rk_steps = 40;
    std::string run_name = "c";
    std::string run_dir = "data";
    int save_each_shifts = 2;  // save confs/measures/seed each
};

#endif  // ECMC_MPI_PARAMS_H
