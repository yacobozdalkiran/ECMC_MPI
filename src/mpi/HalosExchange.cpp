//
// Created by ozdalkiran-l on 1/9/26.
//

#define NDIMS 4

#include "HalosExchange.h"
#include <iostream>

#include "MpiTopology.h"

//Sets the coordinate of the shift in halo. Must be called before the other functions
void mpi::shift::set_coord(Halo &halo, halo_coord coord_) {
    halo.coord = coord_;
}

//Fills the send halo with the corresponding part of the lattice :
//from coord = L - L_shift to coord = L-1 if stype = pos
//from coord = 0 to coord = L_shift-1 if stype = neg
void mpi::shift::fill_halo_send(const GaugeField &field, const GeometryFrozen &geo, Halo &halo,
    shift_type stype) {
    if (halo.coord == X) {
        for (int t = 0; t < geo.L; t++) {
            for (int z = 0; z < geo.L; z++) {
                for (int y = 0; y < geo.L; y++) {
                        for (int x_local = 0; x_local < halo.L_shift; x_local++){
                            int x_glob = x_local;
                            if (stype == pos) x_glob = x_local + geo.L - halo.L_shift;
                            size_t site_glob = geo.index(x_glob, y, z, t);
                            size_t site_loc = halo.index_halo(x_local, y, z, t);
                            for (int mu = 0; mu < 4; mu++) {
                                halo.view_halo_send(site_loc, mu) = field.view_link_const(site_glob, mu);
                            }
                        }
                }
            }
        }
    }

    if (halo.coord == Y) {
        for (int t = 0; t < geo.L; t++) {
            for (int z = 0; z < geo.L; z++) {
                for (int y_local = 0; y_local < halo.L_shift; y_local++) {
                    for (int x = 0; x < geo.L; x++) {
                        int y_glob = y_local;
                        if (stype == pos) y_glob = y_local + geo.L - halo.L_shift;
                        size_t site_glob = geo.index(x, y_glob, z, t);
                        size_t site_loc = halo.index_halo(x, y_local, z, t);
                        for (int mu = 0; mu < 4; mu++) {
                            halo.view_halo_send(site_loc, mu) = field.view_link_const(site_glob, mu);
                        }
                    }
                }
            }
        }
    }
    if (halo.coord == Z) {
        for (int t = 0; t < geo.L; t++) {
            for (int z_local = 0; z_local < halo.L_shift; z_local++) {
                for (int y = 0; y < geo.L; y++) {
                    for (int x = 0; x < geo.L; x++) {
                        int z_glob = z_local;
                        if (stype == pos) z_glob = z_local + geo.L - halo.L_shift;
                        size_t site_glob = geo.index(x, y, z_glob, t);
                        size_t site_loc = halo.index_halo(x, y, z_local, t);
                        for (int mu = 0; mu < 4; mu++) {
                            halo.view_halo_send(site_loc, mu) = field.view_link_const(site_glob, mu);
                        }
                    }
                }
            }
        }
    }
    if (halo.coord == T) {
        for (int t_local= 0; t_local < geo.L; t_local++) {
            for (int z = 0; z < geo.L; z++) {
                for (int y = 0; y < geo.L; y++) {
                    for (int x = 0; x < geo.L; x++) {
                        int t_glob = t_local;
                        if (stype == pos) t_glob = t_local + geo.L - halo.L_shift;
                        size_t site_glob = geo.index(x, y, z, t_glob);
                        size_t site_loc = halo.index_halo(x, y, z, t_local);
                        for (int mu = 0; mu < 4; mu++) {
                            halo.view_halo_send(site_loc, mu) = field.view_link_const(site_glob, mu);
                        }
                    }
                }
            }
        }
    }
}


//Shifts the value of all the links of the lattice of L_shift in direction +coord if stype = pos, -coord if stpye = neg
void mpi::shift::shift_field(GaugeField &field, const mpi::GeometryFrozen &geo, Halo &halo,
    shift_type stype) {
    int L = geo.L;
    int L_shift = halo.L_shift;
    halo_coord coord = halo.coord;
    //We want to avoid cache misses -> we follow the layout of the gauge field and have to write a different loop each time
    if (stype == pos) {
        //The content of the field btw 0 and (L-1-L_shift) -> L_shift - L-1
        if (coord == X) {
            for (int t = 0; t<L; t++) {
                for (int z=  0; z<L; z++) {
                    for (int y=0; y<L; y++) {
                        for (int x = L_shift; x<L; x++) {
                            size_t site = geo.index(x,y,z,t);
                            size_t copied_site = geo.index(x-L_shift, y, z, t);
                            for (int mu = 0; mu<4; mu++) {
                                field.view_link(site, mu) = field.view_link_const(copied_site, mu);
                            }
                        }
                    }
                }
            }
        }

        if (coord == Y) {
            for (int t = 0; t<L; t++) {
                for (int z=  0; z<L; z++) {
                    for (int y=L_shift; y<L; y++) {
                        for (int x = 0; x<L; x++) {
                            size_t site = geo.index(x,y,z,t);
                            size_t copied_site = geo.index(x, y-L_shift, z, t);
                            for (int mu = 0; mu<4; mu++) {
                                field.view_link(site, mu) = field.view_link_const(copied_site, mu);
                            }
                        }
                    }
                }
            }
        }


        if (coord == Z) {
            for (int t = 0; t<L; t++) {
                for (int z= L_shift ; z<L; z++) {
                    for (int y=0; y<L; y++) {
                        for (int x = 0; x<L; x++) {
                            size_t site = geo.index(x,y,z,t);
                            size_t copied_site = geo.index(x, y, z-L_shift, t);
                            for (int mu = 0; mu<4; mu++) {
                                field.view_link(site, mu) = field.view_link_const(copied_site, mu);
                            }
                        }
                    }
                }
            }
        }
        if (coord == T) {
            for (int t = L_shift; t<L; t++) {
                for (int z=  0; z<L; z++) {
                    for (int y=0; y<L; y++) {
                        for (int x = 0; x<L; x++) {
                            size_t site = geo.index(x,y,z,t);
                            size_t copied_site = geo.index(x, y, z, t-L_shift);
                            for (int mu = 0; mu<4; mu++) {
                                field.view_link(site, mu) = field.view_link_const(copied_site, mu);
                            }
                        }
                    }
                }
            }
        }
    }
    if (stype == neg) {
        //The content of the field btw L_shift and (L-1) -> 0 - (L-1-L_shift)
        if (coord == X) {
            for (int t = 0; t<L; t++) {
                for (int z = 0; z<L; z++) {
                    for (int y = 0; y<L; y++) {
                        for (int x = 0; x<L-L_shift; x++) {
                            size_t site = geo.index(x,y,z,t);
                            size_t copied_site = geo.index(x+L_shift,y,z,t);
                            for (int mu = 0; mu<4; mu++) {
                                field.view_link(site, mu) = field.view_link_const(copied_site, mu);
                            }
                        }
                    }
                }
            }
        }

        if (coord == Y) {
            for (int t = 0; t<L; t++) {
                for (int z = 0; z<L; z++) {
                    for (int y = 0; y<L-L_shift; y++) {
                        for (int x = 0; x<L; x++) {
                            size_t site = geo.index(x,y,z,t);
                            size_t copied_site = geo.index(x,y+L_shift,z,t);
                            for (int mu = 0; mu<4; mu++) {
                                field.view_link(site, mu) = field.view_link_const(copied_site, mu);
                            }
                        }
                    }
                }
            }
        }

        if (coord == Z) {
            for (int t = 0; t<L; t++) {
                for (int z = 0; z<L-L_shift; z++) {
                    for (int y = 0; y<L; y++) {
                        for (int x = 0; x<L; x++) {
                            size_t site = geo.index(x,y,z,t);
                            size_t copied_site = geo.index(x,y,z+L_shift,t);
                            for (int mu = 0; mu<4; mu++) {
                                field.view_link(site, mu) = field.view_link_const(copied_site, mu);
                            }
                        }
                    }
                }
            }
        }

        if (coord == T) {
            for (int t = 0; t<L-L_shift; t++) {
                for (int z = 0; z<L; z++) {
                    for (int y = 0; y<L; y++) {
                        for (int x = 0; x<L; x++) {
                            size_t site = geo.index(x,y,z,t);
                            size_t copied_site = geo.index(x,y,z,t+L_shift);
                            for (int mu = 0; mu<4; mu++) {
                                field.view_link(site, mu) = field.view_link_const(copied_site, mu);
                            }
                        }
                    }
                }
            }
        }
    }
}

//Fills the halo_recv of the dest node with the content of the halo_send of the source node
void mpi::shift::exchange_halos(Halo &halo, mpi::MpiTopology &topo, shift_type stype) {
    int source{}, dest{};
    if (halo.coord == X) {
        if (stype == pos) {
            source = topo.x0;
            dest = topo.xL;
        }
        if (stype == neg) {
            source = topo.xL;
            dest = topo.x0;
        }
    }
    if (halo.coord == Y) {
        if (stype == pos) {
            source = topo.y0;
            dest = topo.yL;
        }
        if (stype == neg) {
            source = topo.yL;
            dest = topo.y0;
        }
    }
    if (halo.coord == Z) {
        if (stype == pos) {
            source = topo.z0;
            dest = topo.zL;
        }
        if (stype == neg) {
            source = topo.zL;
            dest = topo.z0;
        }
    }
    if (halo.coord == T) {
        if (stype == pos) {
            source = topo.t0;
            dest = topo.tL;
        }
        if (stype == neg) {
            source = topo.tL;
            dest = topo.t0;
        }
    }
    MPI_Sendrecv(halo.send.data(), 2 * 9 * 4 * halo.V_halo, MPI_DOUBLE, dest, 0,
        halo.recv.data(), 2 * 9 * 4 * halo.V_halo,MPI_DOUBLE, source, 0,
        topo.cart_comm, MPI_STATUS_IGNORE);
}

//Replace the values of the corresponding links of the lattice with those of halo_rec
void mpi::shift::fill_lattice_with_halo_recv(GaugeField &field, const mpi::GeometryFrozen &geo, Halo &halo,
    shift_type stype) {
    int L = geo.L;
    int L_shift = halo.L_shift;
    //If positive shift, the local coordinates of the halo corresponds to the coordinates of the lattice
    if (stype == pos) {
        if (halo.coord == X) {
            for (int t = 0; t<L; t++) {
                for (int z = 0; z<L; z++) {
                    for (int y = 0; y<L; y++) {
                        for (int x = 0; x<L_shift; x++) {
                            size_t index_field = geo.index(x,y,z,t);
                            size_t index_halo = halo.index_halo(x,y,z,t);
                            for (int mu = 0; mu<4; mu++) {
                                field.view_link(index_field, mu) = halo.view_halo_rec_const(index_halo, mu);
                            }
                        }
                    }
                }
            }
        }
        if (halo.coord == Y) {
            for (int t = 0; t<L; t++) {
                for (int z = 0; z<L; z++) {
                    for (int y = 0; y<L_shift; y++) {
                        for (int x = 0; x<L; x++) {
                            size_t index_field = geo.index(x,y,z,t);
                            size_t index_halo = halo.index_halo(x,y,z,t);
                            for (int mu = 0; mu<4; mu++) {
                                field.view_link(index_field, mu) = halo.view_halo_rec_const(index_halo, mu);
                            }
                        }
                    }
                }
            }
        }
        if (halo.coord == Z) {
            for (int t = 0; t<L; t++) {
                for (int z = 0; z<L_shift; z++) {
                    for (int y = 0; y<L; y++) {
                        for (int x = 0; x<L; x++) {
                            size_t index_field = geo.index(x,y,z,t);
                            size_t index_halo = halo.index_halo(x,y,z,t);
                            for (int mu = 0; mu<4; mu++) {
                                field.view_link(index_field, mu) = halo.view_halo_rec_const(index_halo, mu);
                            }
                        }
                    }
                }
            }
        }
        if (halo.coord == T) {
            for (int t = 0; t<L_shift; t++) {
                for (int z = 0; z<L; z++) {
                    for (int y = 0; y<L; y++) {
                        for (int x = 0; x<L; x++) {
                            size_t index_field = geo.index(x,y,z,t);
                            size_t index_halo = halo.index_halo(x,y,z,t);
                            for (int mu = 0; mu<4; mu++) {
                                field.view_link(index_field, mu) = halo.view_halo_rec_const(index_halo, mu);
                            }
                        }
                    }
                }
            }
        }
    }
    //If negative shift, the coordinates of the halo are shifted of L-L_shift to get the correct coordinates on the field
    if (stype == neg) {
        if (halo.coord == X) {
            for (int t = 0; t<L; t++) {
                for (int z = 0; z<L; z++) {
                    for (int y = 0; y<L; y++) {
                        for (int x = 0; x<L_shift; x++) {
                            size_t index_field = geo.index(x+L-L_shift,y,z,t);
                            size_t index_halo = halo.index_halo(x,y,z,t);
                            for (int mu = 0; mu<4; mu++) {
                                field.view_link(index_field, mu) = halo.view_halo_rec_const(index_halo, mu);
                            }
                        }
                    }
                }
            }
        }

        if (halo.coord == Y) {
            for (int t = 0; t<L; t++) {
                for (int z = 0; z<L; z++) {
                    for (int y = 0; y<L_shift; y++) {
                        for (int x = 0; x<L; x++) {
                            size_t index_field = geo.index(x,y+L-L_shift,z,t);
                            size_t index_halo = halo.index_halo(x,y,z,t);
                            for (int mu = 0; mu<4; mu++) {
                                field.view_link(index_field, mu) = halo.view_halo_rec_const(index_halo, mu);
                            }
                        }
                    }
                }
            }
        }

        if (halo.coord == Z) {
            for (int t = 0; t<L; t++) {
                for (int z = 0; z<L_shift; z++) {
                    for (int y = 0; y<L; y++) {
                        for (int x = 0; x<L; x++) {
                            size_t index_field = geo.index(x,y,z+L-L_shift,t);
                            size_t index_halo = halo.index_halo(x,y,z,t);
                            for (int mu = 0; mu<4; mu++) {
                                field.view_link(index_field, mu) = halo.view_halo_rec_const(index_halo, mu);
                            }
                        }
                    }
                }
            }
        }

        if (halo.coord == T) {
            for (int t = 0; t<L_shift; t++) {
                for (int z = 0; z<L; z++) {
                    for (int y = 0; y<L; y++) {
                        for (int x = 0; x<L; x++) {
                            size_t index_field = geo.index(x,y,z,t+L-L_shift);
                            size_t index_halo = halo.index_halo(x,y,z,t);
                            for (int mu = 0; mu<4; mu++) {
                                field.view_link(index_field, mu) = halo.view_halo_rec_const(index_halo, mu);
                            }
                        }
                    }
                }
            }
        }
    }
}

//Performs a shift in every node in direction coord, of type stype (positive of negative) and of length L_shift
void mpi::shift::shift(GaugeField &field, const mpi::GeometryFrozen &geo, Halo &halo, MpiTopology &topo,
    int L_shift, halo_coord coord_, shift_type stype) {
    if (halo.L_shift != L_shift) {
        std::cerr << "Wrong halo size\n";
        return;
    }
    halo.coord = coord_;
    fill_halo_send(field, geo, halo, stype);
    shift_field(field, geo, halo, stype);
    exchange_halos(halo, topo, stype);
    fill_lattice_with_halo_recv(field, geo, halo, stype);
}