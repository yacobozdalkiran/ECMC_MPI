//
// Created by ozdalkiran-l on 1/8/26.
//

#include "GaugeField.h"
#include "../su3/utils.h"


//Initialises the gauge field with random SU3 matrices
void GaugeField::hot_start(std::mt19937_64 &rng) {
    for (size_t site = 0; site < V; site++) {
        for (int mu = 0; mu < 4; mu++) {
            view_link(site, mu) = random_su3(rng);
        }
    }
}

//Initialises the gauge field with identity matrices
void GaugeField::cold_start() {
    for (size_t site = 0; site < V; site++) {
        for (int mu = 0; mu < 4; mu++) {
            view_link(site, mu) = Eigen::Matrix3cd::Identity();
        }
    }
}

//Projects a link on SU3 using Gramm-Schmidt
void GaugeField::projection_su3(size_t site, int mu){
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

