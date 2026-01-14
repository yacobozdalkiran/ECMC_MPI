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

struct RunParams {
    int L_core = 6;
    int n_core_dims = 2;
    bool cold_start = true;
    int L_shift = 2;
    int n_shift = 1; //number of shifts in all the directions
    bool stype_pos=true; //Are the shifts positives ?
    ECMCParams ecmc_params{};
};


#endif //ECMC_MPI_PARAMS_H