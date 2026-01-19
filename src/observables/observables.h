//
// Created by ozdalkiran-l on 1/16/26.
//

#ifndef ECMC_MPI_OBSERVABLES_H
#define ECMC_MPI_OBSERVABLES_H

#include "../gauge/GaugeField.h"

namespace observables {
    double mean_plaquette(const GaugeField &field, const Geometry &geo);
    double mean_plaquette(const GaugeField &field, const mpi::GeometryFrozen &geo);
    double wilson_action(const GaugeField &field, const Geometry &geo);
    void gauge_transform(GaugeField &field, const Geometry &geo, std::mt19937_64 rng);
    double max_drift_det(const GaugeField &field, const Geometry &geo);
}

#endif //ECMC_MPI_OBSERVABLES_H