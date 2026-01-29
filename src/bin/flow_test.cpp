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

int main() {
    //Params
    int L=16;
    int T=32;
    double eps = 0.02;
    int seed = 12345;
    std::mt19937_64 rng(seed);

    //Objects
    Geometry geo(L, T);
    GaugeField field(geo);
    field.hot_start(rng);
    GradientFlow flow(eps, field, geo);

    //Hb steps
    HbParams hp{6.0, 20, 1, 1};
    heatbath::samples(field, geo, hp, rng);

    //Flow
    for (int i = 0; i<1000; i++) {
        std::cout << "t = " << i*eps << ", Q = " << observables::topo_charge_clover(flow.field_c, geo) << "\n";
        flow.rk3_step();
    }
}
