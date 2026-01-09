//
// Created by ozdalkiran-l on 1/9/26.
//

#define NDIMS 4

#include "HalosExchange.h"
#include <iostream>

#include "MpiTopology.h"

//Fills the halo with the cellule defined by mu=0 if i_mu=0 or mu=L-1 if i_mu = 1
void mpi::shift::fill_halo_send(LocalGaugeField &field, const GeometryFrozenMPI &geo, int mu, bool i_mu) {
    if (mu < 0 || mu > 3) std::cerr << "Wrong value of mu\n";
    for (int c1 = 0; c1 < geo.L; c1++) {
        for (int c2 = 0; c2 < geo.L; c2++) {
            for (int c3 = 0; c3 < geo.L; c3++) {
                size_t index_local_halo = geo.index_halo(c1, c2, c3);
                size_t index_global_lattice{};
                if (mu == 0 && i_mu == 0) {
                    index_global_lattice = geo.index(0, c1, c2, c3);
                } else if (mu == 0 && i_mu == 1) {
                    index_global_lattice = geo.index(geo.L - 1, c1, c2, c3);
                } else if (mu == 1 && i_mu == 0) {
                    index_global_lattice = geo.index(c1, 0, c2, c3);
                } else if (mu == 1 && i_mu == 1) {
                    index_global_lattice = geo.index(c1, geo.L - 1, c2, c3);
                } else if (mu == 2 && i_mu == 0) {
                    index_global_lattice = geo.index(c1, c2, 0, c3);
                } else if (mu == 2 && i_mu == 1) {
                    index_global_lattice = geo.index(c1, c2, geo.L - 1, c3);
                } else if (mu == 3 && i_mu == 0) {
                    index_global_lattice = geo.index(c1, c2, c3, 0);
                } else if (mu == 3 && i_mu == 1) {
                    index_global_lattice = geo.index(c1, c2, c3, geo.L - 1);
                }
                for (int nu = 0; nu < NDIMS; nu++) {
                    field.view_halo_send(index_local_halo, nu) = field.view_link_const(index_global_lattice, nu);
                }
            }
        }
    }
}

//Shifts the value of all the links of the lattice in direction +mu
void mpi::shift::shift_pos(mpi::LocalGaugeField &field, const mpi::GeometryFrozenMPI &geo, int mu) {
    int L = geo.L;
    if (mu < 0 || mu > 3) std::cerr << "Wrong value of mu\n";
    for (int cshift = L - 1; cshift > 0; cshift--) {
        for (int c1 = 0; c1 < L; c1++) {
            for (int c2 = 0; c2 < L; c2++) {
                for (int c3 = 0; c3 < L; c3++) {
                    size_t index_site{};
                    size_t index_new{};
                    if (mu == 0) {
                        index_site = geo.index(cshift, c1, c2, c3);
                        index_new = geo.index(cshift - 1, c1, c2, c3);
                    } else if (mu == 1) {
                        index_site = geo.index(c1, cshift, c2, c3);
                        index_new = geo.index(c1, cshift - 1, c2, c3);
                    } else if (mu == 2) {
                        index_site = geo.index(c1, c2, cshift, c3);
                        index_new = geo.index(c1, c2, cshift - 1, c3);
                    } else if (mu == 3) {
                        index_site = geo.index(c1, c2, c3, cshift);
                        index_new = geo.index(c1, c2, c3, cshift - 1);
                    }
                    for (int nu = 0; nu < NDIMS; nu++) {
                        field.view_link(index_site, nu) = field.view_link_const(index_new, nu);
                    }
                }
            }
        }
    }
}

//Shifts the value of all the links of the lattice in direction -mu
void shift_neg(mpi::LocalGaugeField &field, const mpi::GeometryFrozenMPI &geo, int mu) {
    int L = geo.L;
    if (mu < 0 || mu > 3) std::cerr << "Wrong value of mu\n";
    for (int cshift = 0; cshift < L - 1; cshift++) {
        for (int c1 = 0; c1 < L; c1++) {
            for (int c2 = 0; c2 < L; c2++) {
                for (int c3 = 0; c3 < L; c3++) {
                    size_t index_site{};
                    size_t index_new{};
                    if (mu == 0) {
                        index_site = geo.index(cshift, c1, c2, c3);
                        index_new = geo.index(cshift + 1, c1, c2, c3);
                    } else if (mu == 1) {
                        index_site = geo.index(c1, cshift, c2, c3);
                        index_new = geo.index(c1, cshift + 1, c2, c3);
                    } else if (mu == 2) {
                        index_site = geo.index(c1, c2, cshift, c3);
                        index_new = geo.index(c1, c2, cshift + 1, c3);
                    } else if (mu == 3) {
                        index_site = geo.index(c1, c2, c3, cshift);
                        index_new = geo.index(c1, c2, c3, cshift + 1);
                    }
                    for (int nu = 0; nu < NDIMS; nu++) {
                        field.view_link(index_site, nu) = field.view_link_const(index_new, nu);
                    }
                }
            }
        }
    }
}

//Fills the halo_recv of the dest node with the content of the halo_send of the source node
void mpi::shift::exchange_halos(mpi::LocalGaugeField &field, int source, int dest, MPI_Comm comm) {
    MPI_Sendrecv(field.halo_send.data(), 2 * 9 * 4 * field.V_halo, MPI_DOUBLE, dest, 0, field.halo_recv.data(), 2 * 9 * 4 * field.V_halo,
                 MPI_DOUBLE, source, 0, comm, MPI_STATUS_IGNORE);
}

//Replace the values of the corresponding links of the lattice with those of halo_rec
//Halo rec was filled the same way halo_send was, but if halo_send was at coord=L-1, halo_rec contains coord=0
void mpi::shift::fill_lattice_with_halo_recv(mpi::LocalGaugeField &field, const mpi::GeometryFrozenMPI &geo, int mu, bool i_mu) {
    int L = geo.L;
        if (mu < 0 || mu > 3) std::cerr << "Wrong value of mu\n";
        for (int c1 = 0; c1 < L; c1++) {
            for (int c2 = 0; c2 < L; c2++) {
                for (int c3 = 0; c3 < L; c3++) {
                    size_t index_local_halo = geo.index_halo(c1, c2, c3);
                    size_t index_global_lattice{};
                    if (mu == 0 && i_mu == 0) {
                        index_global_lattice = geo.index(L - 1, c1, c2, c3);
                    } else if (mu == 0 && i_mu == 1) {
                        index_global_lattice = geo.index(0, c1, c2, c3);
                    } else if (mu == 1 && i_mu == 0) {
                        index_global_lattice = geo.index(c1, L - 1, c2, c3);
                    } else if (mu == 1 && i_mu == 1) {
                        index_global_lattice = geo.index(c1, 0, c2, c3);
                    } else if (mu == 2 && i_mu == 0) {
                        index_global_lattice = geo.index(c1, c2, L - 1, c3);
                    } else if (mu == 2 && i_mu == 1) {
                        index_global_lattice = geo.index(c1, c2, 0, c3);
                    } else if (mu == 3 && i_mu == 0) {
                        index_global_lattice = geo.index(c1, c2, c3, L - 1);
                    } else if (mu == 3 && i_mu == 1) {
                        index_global_lattice = geo.index(c1, c2, c3, 0);
                    }
                    for (int nu = 0; nu < NDIMS; nu++) {
                        field.view_link(index_global_lattice, nu) = field.view_halo_rec_const(index_local_halo, nu);
                    }
                }
            }
        }
    }

//n shifts in direction mu, c0 et cL are resp. src and dest of exchange_halos
//TODO:check safety for successive shifts -> MPI_Barrier ?
void mpi::shift::n_full_shifts(mpi::LocalGaugeField &field, const mpi::GeometryFrozenMPI &geo, int n, int mu, const mpi::MpiTopology &topo) {
    MPI_Barrier(topo.cart_comm);
    int c0{}, cL{};
    if (mu==0) {
        c0 = topo.x0;
        cL = topo.xL;
    }
    if (mu==1) {
        c0 = topo.y0;
        cL = topo.yL;
    }
    if (mu==2) {
        c0 = topo.z0;
        cL = topo.zL;
    }
    if (mu==3) {
        c0 = topo.t0;
        cL = topo.tL;
    }

    for (int i = 0; i < n; i++) {
        mpi::shift::fill_halo_send(field, geo, mu, true);
        mpi::shift::exchange_halos(field, c0, cL, topo.cart_comm);
        mpi::shift::shift_pos(field, geo, mu);
        mpi::shift::fill_lattice_with_halo_recv(field, geo, mu, true);
    }
    MPI_Barrier(topo.cart_comm);
}
