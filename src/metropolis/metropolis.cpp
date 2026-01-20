//
// Created by ozdalkiran-l on 1/20/26.
//

#include "metropolis.h"
#include <iostream>

#include "../observables/observables.h"
#include "../su3/utils.h"

//Metropolis sweep of the field
void metropolis::sweep(GaugeField &field, const Geometry &geo, const MetroParams &mp,
    size_t &accepted, size_t &proposed, const std::vector<SU3> &set, std::mt19937_64 &rng) {
    SU3 staple;
    SU3 Unew;
    accepted = 0;
    proposed = 0;
    int i_set = 0;
    std::uniform_int_distribution<int> index_set(0, static_cast<int>(set.size()) - 1);
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    for (size_t site = 0; site < geo.V; site++) {
        for (int mu = 0; mu < 4; mu++) {
            field.compute_staple(geo, site, mu, staple);
            auto Umap = field.view_link(site, mu);
            SU3 Uold = Umap;

            for (int i =0; i < mp.N_hits; i++) {
                i_set = index_set(rng);

                Unew = set[i_set] * Uold;

                double old_tr = (Uold * staple).trace().real();
                double new_tr = (Unew * staple).trace().real();
                double dS = -(mp.beta/3.0) * (new_tr - old_tr);

                ++proposed;
                bool accept = false;
                if ((dS <= 0.0)||(unif(rng) < exp(-dS))) accept = true;

                if (accept) {
                    Umap.noalias() = Unew; //No aliasing
                    ++accepted;
                }
            }
        }
    }
}

std::vector<double> metropolis::samples(GaugeField &field, const Geometry &geo, const MetroParams &mp,
    std::mt19937_64 &rng) {

    int size_set = 50;
    auto set = metropolis_set(mp.epsilon,size_set,rng);
    size_t accepted = 0;
    size_t proposed = 0;
    std::cout << "Burn-in...\n";
    for (int i = 0; i < mp.N_burnin; i++) {
        if (i%mp.N_set ==0) set = metropolis_set(mp.epsilon,size_set,rng);
        sweep(field, geo, mp, accepted, proposed, set, rng);
    }
    std::cout << "Finished ! Measuring...\n";
    std::vector<double> measures(mp.N_samples);
    int n_sweeps = mp.N_samples*mp.N_sweeps_meas;
    for (int i = 0; i < n_sweeps; i++) {
        if (i%mp.N_set ==0) set = metropolis_set(mp.epsilon,size_set,rng);
        sweep(field, geo, mp, accepted, proposed, set, rng);
        if (i%mp.N_sweeps_meas==0) {
            auto plaq = observables::mean_plaquette(field, geo);
            measures[i/mp.N_sweeps_meas] = plaq;
            std::cout << "Measure " << i/mp.N_sweeps_meas << ", " << "Metropolis step "<< i <<", <P> = " << plaq << ", acceptance = " << static_cast<double>(accepted) / static_cast<double>(proposed) << std::endl;
            //std::cout << "S = " << wilson_action(links, lat, beta) << endl;
            //double qtopo = topo_charge_clover(links, lat);
            //std::cout << "Qtopo = " << qtopo << endl;
        }
    }
    return measures;
}
