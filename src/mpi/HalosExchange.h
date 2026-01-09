//
// Created by ozdalkiran-l on 1/9/26.
//

#ifndef INC_4D_MPI_HALOSEXCHANGE_H
#define INC_4D_MPI_HALOSEXCHANGE_H

#include "../gauge/LocalGaugeField.h"
#include "../geometry/GeometryFrozen.h"
#include <mpi.h>

namespace mpi::shift {
    void fill_halo_send(LocalGaugeField &field, const GeometryFrozenMPI &geo, int mu, bool i_mu);
    void shift_pos(mpi::LocalGaugeField &field, const mpi::GeometryFrozenMPI &geo, int mu);
    void shift_neg(mpi::LocalGaugeField &field, const mpi::GeometryFrozenMPI &geo, int mu);
    void exchange_halos(mpi::LocalGaugeField &field, int source, int dest, MPI_Comm comm);
    void fill_lattice_with_halo_recv(mpi::LocalGaugeField &field, const mpi::GeometryFrozenMPI &geo, int mu, bool i_mu);
    void n_full_shifts(mpi::LocalGaugeField &field, const mpi::GeometryFrozenMPI &geo, int n, int mu, int c0, int cL, MPI_Comm comm);
};


#endif //INC_4D_MPI_HALOSEXCHANGE_H