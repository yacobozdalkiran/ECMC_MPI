//
// Created by ozdalkiran-l on 1/8/26.
//

#include "GeometryFrozen.h"
#include <cstdint>

mpi::GeometryFrozen::GeometryFrozen(int L_) {
    L = L_;
    V = L*L*L*L;
    neighbors_flat.resize(V*8, SIZE_MAX);
    for (int t = 0; t < L; t++) {
        for (int z = 0; z < L; z++) {
            for (int y = 0; y < L; y++) {
                for (int x = 0; x < L; x++) {
                    size_t i = index(x, y, z, t);
                    if (x+1<L) neighbors_flat[index_neigh(i,0,0)] = index((x + 1), y, z, t);
                    if (x-1 >= 0) neighbors_flat[index_neigh(i,0,1)] = index((x - 1), y, z, t);
                    if (y+1<L) neighbors_flat[index_neigh(i,1,0)] = index(x, (y + 1), z, t);
                    if (y-1 >=0) neighbors_flat[index_neigh(i,1,1)] = index(x, (y - 1), z, t);
                    if (z+1<L) neighbors_flat[index_neigh(i,2,0)] = index(x, y, (z + 1), t);
                    if (z-1 >=0) neighbors_flat[index_neigh(i,2,1)] = index(x, y, (z - 1), t);
                    if (t+1<L) neighbors_flat[index_neigh(i,3,0)] = index(x, y, z, (t + 1));
                    if (t-1>=0) neighbors_flat[index_neigh(i,3,1)] = index(x, y, z, (t - 1));
                }
            }
        }
    }

    links_staples_flat.resize(V*3*6*4, std::make_pair(SIZE_MAX, -1));
    for (int t = 1; t < L-1; t++) {
        for (int z = 1; z < L-1; z++) {
            for (int y = 1; y < L-1; y++) {
                for (int x = 1; x < L-1; x++) {
                    size_t site = index(x, y, z, t); //x
                    for (int mu = 0; mu < 4; mu++) {
                        int j = 0;
                        for (int nu = 0; nu < 4; nu++) {
                            if (nu == mu) continue;

                            size_t xmu = get_neigh(site, mu,0); //x+mu
                            size_t xnu = get_neigh(site, nu,0); //x+nu
                            size_t xmunu = get_neigh(xmu, nu,1); //x+mu-nu
                            size_t xmnu = get_neigh(site, nu,1); //x-nu

                            links_staples_flat[index_staples(site, mu, j, 0)] = {xmu, nu};
                            links_staples_flat[index_staples(site, mu, j, 1)] = {xnu, mu};
                            links_staples_flat[index_staples(site, mu, j, 2)] = {site, nu};

                            links_staples_flat[index_staples(site, mu, j+1, 0)] = {xmunu, nu};
                            links_staples_flat[index_staples(site, mu, j+1, 1)] = {xmnu, mu};
                            links_staples_flat[index_staples(site, mu, j+1, 2)] = {xmnu, nu};

                            j += 2;
                        }
                    }
                }
            }
        }
    }

    is_frozen_flat.resize(V*4, false);

    for (int t = 0; t < L; t++) {
        for (int z = 0; z < L; z++) {
            for (int y = 0; y < L; y++) {
                for (int x = 0; x < L; x++) {
                    size_t site = index(x, y, z, t); //x
                    for (int mu = 0; mu < 4; mu++) {
                        if (x==0||x==L-1||y==0||y==L-1||z==0||z==L-1||t==0||t==L-1) {
                            is_frozen_flat[index_frozen(site, mu)] = true;
                        }
                    }
                }
            }
        }
    }
}
