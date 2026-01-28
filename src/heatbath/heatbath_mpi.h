//
// Created by ozdalkiran-l on 1/28/26.
//

#ifndef ECMC_MPI_HEATBATH_MPI_H
#define ECMC_MPI_HEATBATH_MPI_H

#include "../gauge/GaugeField.h"
#include "../su3/utils.h"
#include "../io/params.h"
#include "../mpi/Halo.h"
#include "../mpi/MpiTopology.h"

namespace mpi::heatbath {
    void hit(GaugeField &field, const GeometryHaloECMC &geo, size_t site, int mu, double beta, SU3 &A, std::mt19937_64 &rng);
    void sweep(GaugeField &field, const GeometryHaloECMC &geo, double beta, int N_hits, std::mt19937_64 &rng);
    std::vector<double> samples(GaugeField &field, const GeometryHaloECMC &geo, HaloObs &obs, MpiTopology &topo, const HbParams &params, std::mt19937_64 &rng);
}

#endif //ECMC_MPI_HEATBATH_MPI_H