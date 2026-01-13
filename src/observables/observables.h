//
// Created by ozdalkiran-l on 1/8/26.
//

#ifndef INC_4D_MPI_OBSERVABLES_H
#define INC_4D_MPI_OBSERVABLES_H

#include "../gauge/GaugeField.h"
#include "../mpi/Halo.h"
#include "../mpi/MpiTopology.h"


namespace observables {
    double mean_plaquette(const GaugeField &field, const Geometry &geo);
    double mean_plaquette(const GaugeField &field, const mpi::GeometryFrozen &geo);
    double wilson_action(const GaugeField &field, const Geometry &geo);
}

namespace mpi::observables {
    double sum_plaquette_bulk(const GaugeField &field, const mpi::GeometryFrozen &geo);
    double sum_plaquette_boundaries(const GaugeField &field, const mpi::GeometryFrozen &geo, const HaloObs &halo_obs);
    double mean_plaquette_local(const GaugeField &field, const mpi::GeometryFrozen &geo, HaloObs &halo_obs, mpi::MpiTopology &topo);
    double mean_plaquette_global(const GaugeField &field, const mpi::GeometryFrozen &geo, HaloObs &halo_obs, mpi::MpiTopology &topo);
}

#endif //INC_4D_MPI_OBSERVABLES_H