//
// Created by ozdalkiran-l on 1/8/26.
//

#include "utils.h"

#include <iostream>
#include "../gauge/GaugeField.h"

//Generates a random SU3 matrix using QR decomposition
SU3 random_su3(std::mt19937_64 &rng) {
    std::normal_distribution<double> gauss(0.0, 1.0);
    Eigen::Matrix3cd z;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            z(i,j) = Complex(gauss(rng), gauss(rng));

    //QR Decomposition
    Eigen::HouseholderQR<SU3> qr(z);
    SU3 Q = qr.householderQ();

    Complex detQ = Q.determinant();
    Q /= std::pow(detQ, 1.0/3.0);

    return Q;
}

//Embbeds a SU2 matrix into a SU3 matrix according to subgroup (i,j)
SU3 su2_quaternion_to_su3(const std::array<double,4> &su2, int i, int j){
    if (i==j) std::cerr<<"i = j wrong embedding\n";
    SU3 X;
    int k = 3-i-j;
    X.setZero();
    X(k,k) = Complex(1.0,0.0);
    X(i,i) = Complex(su2[0],su2[3]);
    X(j,j) = Complex(su2[0],-su2[3]);
    X(i,j) = Complex(su2[2],su2[1]);
    X(j,i) = Complex(-su2[2],su2[1]);
    return X;
}

//Returns a random SU3 matrix epsilon close to identity
SU3 random_SU3_epsilon(double epsilon, std::mt19937_64 &rng) {
    std::uniform_real_distribution<double> unif(-0.5,0.5);
    std::array<double,4> x = {0.0, 0.0, 0.0, 0.0};
    SU3 M = SU3::Identity();

    //double r0 = unif(rng);
    double r1 = unif(rng);
    double r2 = unif(rng);
    double r3 = unif(rng);
    double norm = sqrt(r1*r1 + r2*r2 + r3*r3);
    x[0] = sqrt(1-epsilon*epsilon);
    x[1] = epsilon * r1 / norm;
    x[2] = epsilon * r2 / norm;
    x[3] = epsilon * r3 / norm;
    M *= su2_quaternion_to_su3(x, 0,1);

    //r0 = unif(rng);
    r1 = unif(rng);
    r2 = unif(rng);
    r3 = unif(rng);
    norm = sqrt(r1*r1 + r2*r2 + r3*r3);
    x[0] = sqrt(1-epsilon*epsilon);
    x[1] = epsilon * r1 / norm;
    x[2] = epsilon * r2 / norm;
    x[3] = epsilon * r3 / norm;
    M *= su2_quaternion_to_su3(x, 0,2);

    //r0 = unif(rng);
    r1 = unif(rng);
    r2 = unif(rng);
    r3 = unif(rng);
    norm = sqrt(r1*r1 + r2*r2 + r3*r3);
    x[0] = sqrt(1-epsilon*epsilon);
    x[1] = epsilon * r1 / norm;
    x[2] = epsilon * r2 / norm;
    x[3] = epsilon * r3 / norm;
    M *= su2_quaternion_to_su3(x, 1,2);

    return M;
}

//Creates a set of SU3 matrices close to identity
std::vector<SU3> ecmc_set(double epsilon, std::vector<SU3> &set, std::mt19937_64 &rng) {
    size_t size = set.size()-1;
    set[0] = SU3::Identity();
    for (size_t i = 1; i < size+1; i+=2) {
        set[i] = random_SU3_epsilon(epsilon, rng);
        set[i+1] = set[i].adjoint();
    }
    return set;
}


