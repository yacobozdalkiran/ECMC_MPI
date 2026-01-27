//
// Created by ozdalkiran-l on 1/16/26.
//

#ifndef ECMC_MPI_OBSERVABLES_H
#define ECMC_MPI_OBSERVABLES_H

#include "../gauge/GaugeField.h"

namespace observables {
    double mean_plaquette(const GaugeField &field, const Geometry &geo);
    double wilson_action(const GaugeField &field, const Geometry &geo);
    void gauge_transform(GaugeField &field, const Geometry &geo, std::mt19937_64 rng);
    double max_drift_det(const GaugeField &field, const Geometry &geo);
    SU3 clover_site(const GaugeField &field, const Geometry &geo, size_t site, int mu, int nu);
    inline int levi_civita(int mu, int nu, int rho, int sigma) {
        if (mu == nu || mu == rho || mu == sigma || nu == rho || nu == sigma || rho == sigma) return 0;
        int inv = 0;
        inv += (mu > nu);
        inv += (mu > rho);
        inv += (mu > sigma);
        inv += (nu > rho);
        inv += (nu > sigma);
        inv += (rho > sigma);

        return (inv & 1) ? -1 : +1; //Parity test on inv
    }
    double local_topo_charge_clover(const GaugeField &field, const Geometry &geo, size_t site);
    double topo_charge_clover(const GaugeField &field, const Geometry &geo);
}

#endif //ECMC_MPI_OBSERVABLES_H