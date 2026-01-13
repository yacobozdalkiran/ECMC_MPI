//
// Created by ozdalkiran-l on 1/12/26.
//

#ifndef INC_4D_MPI_HALO_H
#define INC_4D_MPI_HALO_H

#include "../gauge/GaugeField.h"

enum halo_coord {
    X,
    Y,
    Z,
    T
};


class Halo {
public:
    std::vector<Complex> send;
    std::vector<Complex> recv;
    halo_coord coord; //Direction of the shift
    int L_shift; //Length of the halo, such that V_halo = L*L*L*L_shift
    int L; //Size of the square lattice
    int V_halo; //Volume (number of sites) of the halos

    //Creates a shift halo of size size_ for a square gauge configuration of size L_**4
    //halo_coord is the axis of the shift (X,Y,Z,T)
    explicit Halo(int L_shift_, int L_, halo_coord coord_) {
        L_shift = L_shift_;
        L = L_;
        coord = coord_;
        V_halo = L*L*L*L_shift;
        send.resize(V_halo*4*9);
        recv.resize(V_halo*4*9);
    }

    //Index function for local coordinates x,y,z,t of halo send/recv
    [[nodiscard]] size_t index_halo(int x, int y, int z, int t) const {
        size_t index{};
        if (coord == X)
            index = ((static_cast<size_t>(t)*L + z) * L + y)*L_shift+ x;
        if (coord == Y)
            index = ((static_cast<size_t>(t)*L + z) * L_shift+ y)* L + x;
        if (coord == Z)
            index =  ((static_cast<size_t>(t)*L_shift+ z) * L+ y)* L + x;
        if (coord == T)
            index =  ((static_cast<size_t>(t)*L + z) * L + y)* L + x;
        return index;
    }

    //Non const mapping of halo_send to SU3 matrices
    Eigen::Map<SU3> view_halo_send(size_t site, int mu) {
        return Eigen::Map<SU3>(&send[(site * 4 + mu) * 9]);
    }

    //Const mapping of halo_send to SU3 matrices
    [[nodiscard]] Eigen::Map<const SU3> view_halo_send_const(size_t site, int mu) const {
        return Eigen::Map<const SU3>(&send[(site * 4 + mu) * 9]);
    }

    //Non const mapping of halo_recv to SU3 matrices
    Eigen::Map<SU3> view_halo_rec(size_t site, int mu) {
        return Eigen::Map<SU3>(&recv[(site * 4 + mu) * 9]);
    }

    //Const mapping of halo_recv to SU3 matrices
    [[nodiscard]] Eigen::Map<const SU3> view_halo_rec_const(size_t site, int mu) const {
        return Eigen::Map<const SU3>(&recv[(site * 4 + mu) * 9]);
    }
};

#endif //INC_4D_MPI_HALO_H