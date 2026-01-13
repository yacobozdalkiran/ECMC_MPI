//
// Created by ozdalkiran-l on 1/9/26.
//

#ifndef INC_4D_MPI_HALOSEXCHANGE_H
#define INC_4D_MPI_HALOSEXCHANGE_H

#include "../geometry/GeometryFrozen.h"
#include "MpiTopology.h"
#include "Halo.h"

enum shift_type {
    pos, //if the shift is in direction +coord
    neg //if the shift is in direction -coord
};

namespace mpi::shift {
    void set_coord(Halo &halo, halo_coord coord_);
    void fill_halo_send(const GaugeField &field, const GeometryFrozen &geo, Halo &halo, shift_type stype);
    void shift_field(GaugeField &field, const mpi::GeometryFrozen &geo, Halo &halo, shift_type stype);
    void exchange_halos(Halo &halo, mpi::MpiTopology &topo, shift_type stype, MPI_Request* req);
    void fill_lattice_with_halo_recv(GaugeField &field, const mpi::GeometryFrozen &geo, Halo& halo, shift_type stype);
    void shift(GaugeField &field, const mpi::GeometryFrozen &geo, Halo &halo, MpiTopology &topo, int L_shift, halo_coord coord_, shift_type stype);
};

namespace mpi::observables {
    void fill_halo_obs_send(const GaugeField &field, const mpi::GeometryFrozen &geo, HaloObs &halo_obs);
    void exchange_halos_obs(HaloObs &halo_obs, mpi::MpiTopology &topo, MPI_Request* reqs);
}


#endif //INC_4D_MPI_HALOSEXCHANGE_H