//
// Created by ozdalkiran-l on 1/8/26.
//

#ifndef INC_4D_MPI_GAUGEFIELD_H
#define INC_4D_MPI_GAUGEFIELD_H

#include <Eigen/Dense>
#include <vector>
#include <complex>
#include <random>
#include "../geometry/Geometry.h"

using Complex = std::complex<double>;
using SU3 = Eigen::Matrix3cd;

class GaugeField {
    int L;
    int T;
    size_t V;
    std::vector<Complex> links;

public:

    //Initialises a cold gauge conf
    explicit GaugeField(int L_, int T_):L(L_), T(T_), V(L_*L_*L_*T_), links(V*4*9, Complex(0.0,0.0)) {
        for (size_t site = 0; site < V; site++) {
            for (int mu = 0; mu<4; mu++) {
                view_link(site, mu) = SU3::Identity();
            }
        }
    }

    void hot_start(std::mt19937_64 &rng);
    void cold_start();

    //Non const mapping of links to SU3 matrices
    Eigen::Map<SU3> view_link(size_t site, int mu) {
        return Eigen::Map<SU3>(&links[(site * 4 + mu) * 9]);
    }
    //Const mapping of links to SU3 matrices
    Eigen::Map<const SU3> view_link_const(size_t site, int mu) const{
        return Eigen::Map<const SU3>(&links[(site * 4 + mu) * 9]);
    }

    void projection_su3(size_t site, int mu);
};

void compute_staple(const GaugeField &field, const Geometry &geo, size_t site, int mu, SU3 &staple);

#endif //INC_4D_MPI_GAUGEFIELD_H