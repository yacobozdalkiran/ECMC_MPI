//
// Created by ozdalkiran-l on 1/8/26.
//

#ifndef INC_4D_MPI_ECMCUPDATER_H
#define INC_4D_MPI_ECMCUPDATER_H

#include <random>
#include "../gauge/GaugeField.h"


class EcmcUpdater {
private:
    std::mt19937_64 rng;
    std::uniform_real_distribution<double> unif01;

public:
    //Initializes the random generators with random seed
    EcmcUpdater(): rng(std::random_device{}()), unif01(0.0,1.0) {
    }
    //Initializes the random generators with specific seed
    EcmcUpdater(size_t seed): rng(seed), unif01(0.0,1.0) {
    }

    static void compute_list_staples(const GaugeField &field, const Geometry &geo, size_t site, int mu, std::array<SU3,6> &list_staple);
    static void solve_reject(double A, double B, double &gamma, double &reject, int epsilon);
    void compute_reject_angles(const GaugeField &field, size_t site, int mu, const std::array<SU3,6> &list_staple,
        const SU3 &R, int epsilon, const double &beta, std::array<double,6> &reject_angles);
    size_t selectVariable(const std::array<double,4> &probas);
    std::pair<std::pair<size_t, int>,int> lift_improved(const GaugeField &field, const Geometry &geo, size_t site, int mu, int j,
        SU3 &R, const SU3 &lambda_3, const std::vector<SU3> &set);
    static void ecmc_update(GaugeField &field, size_t site, int mu, double theta, int epsilon, const SU3 &R);
    std::vector<double> ecmc_samples_improved(GaugeField &field, const Geometry &geo, double beta, int N_samples,
    double param_theta_sample, double param_theta_refresh, bool poisson, double epsilon_set);
};

//Returns the sign of the double x
inline int dsign(double x) {
    //fonction signe pour double
    if (x>0) return 1;
    if (x<0) return -1;
    return 0;
}

#endif //INC_4D_MPI_ECMCUPDATER_H