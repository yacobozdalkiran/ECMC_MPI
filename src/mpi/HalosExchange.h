//
// Created by ozdalkiran-l on 1/9/26.
//

#ifndef INC_4D_MPI_HALOSEXCHANGE_H
#define INC_4D_MPI_HALOSEXCHANGE_H

#include "../geometry/GeometryFrozen.h"
#include "MpiTopology.h"
#include "Halo.h"
#include "../io/params.h"


namespace mpi::shift {
    void fill_halo_send(const GaugeField &field, const GeometryFrozen &geo, Halo &halo, const ShiftParams &sp);
    void shift_field(GaugeField &field, const mpi::GeometryFrozen &geo, Halo &halo, const ShiftParams &sp);
    void exchange_halos(Halo &halo, mpi::MpiTopology &topo, const ShiftParams &sp, MPI_Request* req);
    void fill_lattice_with_halo_recv(GaugeField &field, const mpi::GeometryFrozen &geo, Halo& halo, const ShiftParams &sp);
    void shift(GaugeField &field, const mpi::GeometryFrozen &geo, Halo &halo, MpiTopology &topo, const ShiftParams &sp);
};

namespace mpi::observables {
    void fill_halo_obs_send(const GaugeField &field, const mpi::GeometryFrozen &geo, HaloObs &halo_obs);
    void exchange_halos_obs(HaloObs &halo_obs, mpi::MpiTopology &topo, MPI_Request* reqs);
}

namespace mpi::ecmc {
    void fill_halos_ecmc(const GaugeField &field, const mpi::GeometryFrozen &geo, HaloECMC &halo);
    void exchange_halos_ecmc(const mpi::GeometryFrozen &geo, HaloECMC &halo, mpi::MpiTopology &topo);
}

#endif //INC_4D_MPI_HALOSEXCHANGE_H