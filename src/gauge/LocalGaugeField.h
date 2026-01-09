//
// Created by ozdalkiran-l on 1/9/26.
//

#ifndef INC_4D_MPI_LOCALGAUGEFIELD_H
#define INC_4D_MPI_LOCALGAUGEFIELD_H

#include <complex>
#include <Eigen/Dense>
#include <random>
#include "../geometry/GeometryFrozen.h"

using Complex = std::complex<double>;
using SU3 = Eigen::Matrix3cd;

namespace mpi {
    class LocalGaugeField {
    public:
        int L;
        size_t V;
        size_t V_halo;
        std::vector<Complex> links; //Links of the gauge field
        std::vector<Complex> halo_send; //Halo to send links
        std::vector<Complex> halo_recv; //Halo to receive links

        //Initialises a cold gauge conf and the halos to zeros
        explicit LocalGaugeField(int L_):
        L(L_), V(L_*L_*L_*L_), V_halo(L_*L_*L_),
        links(V*4*9, Complex(0.0,0.0)),
        halo_send(V_halo*4*9, Complex(0.0, 0.0)),
        halo_recv(V_halo*4*9, Complex(0.0, 0.0))
        {
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
        [[nodiscard]] Eigen::Map<const SU3> view_link_const(size_t site, int mu) const{
            return Eigen::Map<const SU3>(&links[(site * 4 + mu) * 9]);
        }

        void projection_su3(size_t site, int mu);

        //Non const mapping of halo_send to SU3 matrices
        Eigen::Map<SU3> view_halo_send(size_t site, int mu) {
            return Eigen::Map<SU3>(&halo_send[(site * 4 + mu) * 9]);
        }

        //Const mapping of halo_send to SU3 matrices
        [[nodiscard]] Eigen::Map<const SU3> view_halo_send_const(size_t site, int mu) const {
            return Eigen::Map<const SU3>(&halo_send[(site * 4 + mu) * 9]);
        }

        //Non const mapping of halo_recv to SU3 matrices
        Eigen::Map<SU3> view_halo_rec(size_t site, int mu) {
            return Eigen::Map<SU3>(&halo_recv[(site * 4 + mu) * 9]);
        }

        //Const mapping of halo_recv to SU3 matrices
        [[nodiscard]] Eigen::Map<const SU3> view_halo_rec_const(size_t site, int mu) const {
            return Eigen::Map<const SU3>(&halo_recv[(site * 4 + mu) * 9]);
        }
    };

    namespace ecmc {
        void compute_staple(const mpi::LocalGaugeField &field, const mpi::GeometryFrozenMPI &geo, size_t site, int mu, SU3 &staple);
    }
}

#endif //INC_4D_MPI_LOCALGAUGEFIELD_H