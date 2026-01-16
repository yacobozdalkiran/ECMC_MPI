//
// Created by ozdalkiran-l on 1/8/26.
//

#ifndef INC_4D_MPI_ECMCUPDATER_H
#define INC_4D_MPI_ECMCUPDATER_H

#include <random>
#include "../gauge/GaugeField.h"
#include "../io/params.h"

namespace ecmc {
    void compute_list_staples(const GaugeField &field, const Geometry &geo, size_t site, int mu, std::array<SU3,6> &list_staple);
    void solve_reject(double A, double B, double &gamma, double &reject, int epsilon);
    void compute_reject_angles(const GaugeField &field, size_t site, int mu, const std::array<SU3,6> &list_staple,
        const SU3 &R, int epsilon, const double &beta, std::array<double,6> &reject_angles, std::mt19937_64 &rng);
    size_t selectVariable(const std::array<double,4> &probas, std::mt19937_64 &rng);
    std::pair<std::pair<size_t, int>,int> lift_improved(const GaugeField &field, const Geometry &geo, size_t site, int mu, int j,
        SU3 &R, const SU3 &lambda_3, const std::vector<SU3> &set, std::mt19937_64 &rng);
    void update(GaugeField &field, size_t site, int mu, double theta, int epsilon, const SU3 &R);
    std::vector<double> samples_improved(GaugeField &field, const Geometry &geo, const ECMCParams &params, std::mt19937_64 &rng);
}

#endif //INC_4D_MPI_ECMCUPDATER_H