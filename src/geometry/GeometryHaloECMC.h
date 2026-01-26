//
// Created by ozdalkiran-l on 1/26/26.
//

#ifndef ECMC_MPI_GEOMETRYHALOECMC_H
#define ECMC_MPI_GEOMETRYHALOECMC_H

#include <vector>
#include <array>
#include <iostream>

enum dir {
    up,
    down
};

class GeometryHaloECMC {

public:
    int L;
    size_t V;
    size_t V_halo;
private:
    std::vector<size_t> neighbors;

public:
    explicit GeometryHaloECMC(int L_);

    //Index function for links
    [[nodiscard]] size_t index(int x, int y, int z, int t) const {
        return ((static_cast<size_t>(t)*L + z) * L + y)*L + x;
    }

    //Returns the index of the halo site corresponding to (c1,c2,c3) with order x,y,z,t
    [[nodiscard]] size_t index_halo_ecmc(int c1, int c2, int c3) const {
        return (c3*L + c2)*L+c1;
    }

    [[nodiscard]] size_t index_site_neigh(int x, int y, int z, int t) const {
        if (x==L) {
            return (V+index_halo_ecmc(y,z,t));
        }
        if (y==L) {
            return (V+V_halo+index_halo_ecmc(x,z,t));
        }
        if (z==L) {
            return (V+2*V_halo+index_halo_ecmc(x,y,t));
        }
        if (t==L) {
            return (V+3*V_halo+index_halo_ecmc(x,y,z));
        }
        return index(x,y,z,t);
    }

    [[nodiscard]] size_t index_neigh(size_t site, int mu, dir d) {
        return site*8+mu*2+d;
    }
};


#endif //ECMC_MPI_GEOMETRYHALOECMC_H