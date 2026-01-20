//
// Created by ozdalkiran-l on 1/20/26.
//

#ifndef ECMC_MPI_METROPOLIS_H
#define ECMC_MPI_METROPOLIS_H

#include "../gauge/GaugeField.h"
#include "../io/params.h"

namespace metropolis {
    void sweep(GaugeField &field, const Geometry &geo, const MetroParams &mp, size_t &accepted, size_t &proposed, const std::vector<SU3> &set, std::mt19937_64 &rng);
    std::vector<double> samples(GaugeField &field, const Geometry &geo, const MetroParams &mp, std::mt19937_64 &rng);
}

#endif //ECMC_MPI_METROPOLIS_H