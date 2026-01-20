//
// Created by ozdalkiran-l on 1/14/26.
//

#ifndef ECMC_MPI_PARAMS_H
#define ECMC_MPI_PARAMS_H


struct ECMCParams {
    double beta = 6.0;
    int N_samples = 10;
    double param_theta_sample = 100;
    double param_theta_refresh = 30;
    bool poisson = false;
    double epsilon_set = 0.15;
};

struct MetroParams {
    double beta = 6.0;
    double epsilon = 0.15;
    int N_set = 50; //Number of sweeps between set change
    int N_samples = 10; //Number of measures (samples)
    int N_sweeps_meas = 1; //Number of sweeps between measures
    int N_hits = 1;
    int N_burnin = 0;
};

struct RunParams {
    int L_core = 6;
    int n_core_dims = 2;
    bool cold_start = true;
    int L_shift = 2;
    int n_shift = 1; //number of shifts in all the directions
    bool stype_pos=true; //Are the shifts positives ?
    ECMCParams ecmc_params{};
};

struct RunParamsSC {
    int L=4;
    int T=4;
    bool cold_start=true;
    ECMCParams ecmc_params{};
};

struct RunParamsMetro {
    int L = 4;
    int T = 4;
    bool cold_start=true;
    MetroParams mp{};
};

#endif //ECMC_MPI_PARAMS_H