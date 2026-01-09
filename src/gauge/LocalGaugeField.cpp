//
// Created by ozdalkiran-l on 1/9/26.
//

#include "LocalGaugeField.h"
#include "../su3/utils.h"


//Initialises the gauge field with random SU3 matrices
void mpi::LocalGaugeField::hot_start(std::mt19937_64 &rng) {
    for (size_t site = 0; site < V; site++) {
        for (int mu = 0; mu < 4; mu++) {
            view_link(site, mu) = random_su3(rng);
        }
    }
}

//Initialises the gauge field with identity matrices
void mpi::LocalGaugeField::cold_start() {
    for (size_t site = 0; site < V; site++) {
        for (int mu = 0; mu < 4; mu++) {
            view_link(site, mu) = Eigen::Matrix3cd::Identity();
        }
    }
}

//Projects a link on SU3 using Gramm-Schmidt
void mpi::LocalGaugeField::projection_su3(size_t site, int mu){
    SU3 U = view_link(site, mu);
    Eigen::Vector3cd c0 = U.col(0);
    Eigen::Vector3cd c1 = U.col(1);
    Eigen::Vector3cd c2 = U.col(2);

    c0.normalize();
    c1 -= c0 * (c0.adjoint() * c1)(0,0);
    c1.normalize();
    c2 -= c0 * (c0.adjoint() * c2)(0,0);
    c2 -= c1 * (c1.adjoint() * c2)(0,0);
    c2.normalize();

    U.col(0) = c0;
    U.col(1) = c1;
    U.col(2) = c2;

    // Fix determinant to 1
    Complex det = U.determinant();
    U /= exp(Complex(0.0, arg(det) / 3.0));
}

//Computes the sum of staples of link (site, mu)
void mpi::ecmc::compute_staple(const mpi::LocalGaugeField &field, const mpi::GeometryFrozenMPI &geo, size_t site, int mu, SU3 &staple) {
    staple.setZero();
    int j = 0;
    while (j < 6) {

        auto U0 = field.view_link_const(geo.get_link_staple(site,mu,j,0).first,geo.get_link_staple(site,mu,j,0).second);
        auto U1 = field.view_link_const(geo.get_link_staple(site,mu,j,1).first, geo.get_link_staple(site,mu,j,1).second);
        auto U2 = field.view_link_const(geo.get_link_staple(site,mu,j,2).first, geo.get_link_staple(site,mu,j,2).second);
        staple += U0 * U1.adjoint() * U2.adjoint();

        auto V0 = field.view_link_const(geo.get_link_staple(site,mu,j+1,0).first, geo.get_link_staple(site,mu,j+1,0).second);
        auto V1 = field.view_link_const(geo.get_link_staple(site,mu,j+1,1).first, geo.get_link_staple(site,mu,j+1,1).second);
        auto V2 = field.view_link_const(geo.get_link_staple(site,mu,j+1,2).first, geo.get_link_staple(site,mu,j+1,2).second);
        staple += V0.adjoint() * V1.adjoint() * V2;
        j += 2;
    }
}
