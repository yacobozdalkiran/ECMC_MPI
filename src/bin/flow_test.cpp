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
    double tlim = 3.0;
    for (int i = 0; i<(tlim+0.02)/0.02; i++) {
        flow.rk3_step();
    }
    auto qe = observables::topo_q_e_clover(flow.field_c, geo);
    std::cout
    << "t = " << io::format_double(tlim, precision_t)
    << ", Q = " << io::format_double(qe.first, precision_Q)
    << ", t^2*E = " << io::format_double(pow(tlim,2)*qe.second, precision_tsqE)
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

    //ECMC params and thermalisation
    std::cout << "Thermalisation...\n";
    HbParams p{6.0, 100, 3, 1};
    heatbath::samples(field, geo, p, rng);

    //Samples and Flow
    int N = 5;
    p.N_samples = 20; //tau_int = 2 for the plaquette
    for (int count = 0; count<N; count++) {
        heatbath::samples(field, geo, p, rng);
        if (count==4) flow_charge(field, geo, eps);
    }
    std::string filename = "Q0.9";
    io::ildg::save_ildg(field, geo, filename);
    std::cout << field.view_link_const(0,0) << "\n";
}
