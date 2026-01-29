//
// Created by ozdalkiran-l on 1/29/26.
//

#include <iostream>
#include "../gauge/GaugeField.h"
#include "../geometry/Geometry.h"
#include "../flow/gradient_flow.h"
#include "../heatbath/heatbath.h"
#include "../observables/observables.h"
#include "../io/params.h"
#include "../io/io.h"

int main() {
    //Params
    int L=8;
    int T=16;
    double eps = 0.02;
    int seed = 12345;
    std::mt19937_64 rng(seed);

    //Objects
    Geometry geo(L, T);
    GaugeField field(geo);
    field.hot_start(rng);

    //Hb steps
    HbParams hp{6.0, 50, 1, 1};
    heatbath::samples(field, geo, hp, rng);

    //Flow
    int precision_Q = 6;
    int precision_tsqE = 3;
    int precision_t = 2;

    GradientFlow flow(eps, field, geo);
    for (int i = 0; i<1000; i++) {
        auto qe = observables::topo_q_e_clover(flow.field_c, geo);
        std::cout
        << "t = " << io::format_double(i*eps, precision_t)
        << ", Q = " << io::format_double(qe.first, precision_Q)
        << ", t^2*E = " << io::format_double(eps*i*eps*i*qe.second, precision_tsqE)
        << ", P = " << io::format_double(observables::mean_plaquette(flow.field_c, geo), precision_Q)
        << "\n";
        flow.rk3_step();
        auto M = flow.field_c.view_link(0,0);
        std::cout << (M * M.adjoint()).trace() << std::endl;
    }
}
