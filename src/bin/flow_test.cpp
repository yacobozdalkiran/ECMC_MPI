//
// Created by ozdalkiran-l on 1/29/26.
//

#include <iostream>

#include "../ecmc/ecmc.h"
#include "../gauge/GaugeField.h"
#include "../geometry/Geometry.h"
#include "../flow/gradient_flow.h"
#include "../heatbath/heatbath.h"
#include "../observables/observables.h"
#include "../io/params.h"
#include "../io/io.h"

void flow_charge(const GaugeField &field, const Geometry &geo, double eps) {
    //Flow
    int precision_Q = 6;
    int precision_tsqE = 3;
    int precision_t = 2;
    GradientFlow flow(eps, field, geo);
    for (int i = 0; i<(2.5/0.02); i++) {
        flow.rk3_step();
    }
    auto qe = observables::topo_q_e_clover(flow.field_c, geo);
    std::cout
    << "t = " << io::format_double(2.48, precision_t)
    << ", Q = " << io::format_double(qe.first, precision_Q)
    << ", t^2*E = " << io::format_double(2.48*2.48*qe.second, precision_tsqE)
    << ", P = " << io::format_double(observables::mean_plaquette(flow.field_c, geo), precision_Q)
    << "\n";
}

int main() {
    //Params
    int L=8;
    int T=8;
    double eps = 0.02;
    int seed = 12345;
    std::mt19937_64 rng(seed);

    //Objects
    Geometry geo(L, T);
    GaugeField field(geo);
    field.hot_start(rng);

    //Hb params and thermalisation
    std::cout << "Thermalisation...\n";
    ECMCParams ep{6.0, 500, 6050, 1000, false, 0.15};
    ecmc::samples_improved(field, geo, ep, rng);

    //Samples and Flow
    int N = 10;
    ep.N_samples = 100; //tau_int = 20 for the plaquette
    for (int count = 0; count<N; count++) {
        ecmc::samples_improved(field, geo, ep, rng);
        flow_charge(field, geo, eps);
    }
}
