//
// Created by ozdalkiran-l on 1/29/26.
//

#include <iostream>
#include "../su3/utils.h"


int main() {
    int seed = 1234;
    std::mt19937_64 rng(seed);
    SU3 M = random_su3(rng);
    std::cout << M << "\n";
    proj_lie_su3(M);
    std::cout << "M :\n" << M <<"\n";
    std::cout << "M_adj :\n" << M.adjoint() << ",\n Trace(M) : " << M.trace() << "\n";
    SU3 expM = exp_analytic(M, 1.0);
    std::cout << expM << "\n";
    std::cout << "Determinant : " << expM.determinant() << ", exp(M)*adj :\n" << expM.adjoint()*expM << "\n";
}