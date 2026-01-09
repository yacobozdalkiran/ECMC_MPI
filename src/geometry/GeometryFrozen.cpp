//
// Created by ozdalkiran-l on 1/8/26.
//

#include "GeometryFrozen.h"

GeometryFrozen::GeometryFrozen(int L_, int T_) {
    L = L_;
    T = T_;
    V = T*L*L*L;
    neighbors_flat.resize(V*8);
    for (int t = 0; t < T; t++) {
        for (int z = 0; z < L; z++) {
            for (int y = 0; y < L; y++) {
                for (int x = 0; x < L; x++) {
                    size_t i = index(x, y, z, t);
                    neighbors_flat[index_neigh(i,0,0)] = index((x + 1) % L, y, z, t);
                    neighbors_flat[index_neigh(i,0,1)] = index((x - 1 + L) % L, y, z, t);
                    neighbors_flat[index_neigh(i,1,0)] = index(x, (y + 1) % L, z, t);
                    neighbors_flat[index_neigh(i,1,1)] = index(x, (y - 1 + L) % L, z, t);
                    neighbors_flat[index_neigh(i,2,0)] = index(x, y, (z + 1) % L, t);
                    neighbors_flat[index_neigh(i,2,1)] = index(x, y, (z - 1 + L) % L, t);
                    neighbors_flat[index_neigh(i,3,0)] = index(x, y, z, (t + 1) % T);
                    neighbors_flat[index_neigh(i,3,1)] = index(x, y, z, (t - 1 + T) % T);
                }
            }
        }
    }

    links_staples_flat.resize(V*3*6*4);
    for (int t = 0; t < T; t++) {
        for (int z = 0; z < L; z++) {
            for (int y = 0; y < L; y++) {
                for (int x = 0; x < L; x++) {
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

    for (int t = 0; t < T; t++) {
        for (int z = 0; z < L; z++) {
            for (int y = 0; y < L; y++) {
                for (int x = 0; x < L; x++) {
                    size_t site = index(x, y, z, t); //x
                    for (int mu = 0; mu < 4; mu++) {
                        if (x==0||x==L-1||y==0||y==L-1||z==0||z==L-1||t==0||t==T-1) {
                            is_frozen_flat[index_frozen(site, mu)] = true;
                        }
                    }
                }
            }
        }
    }
}

mpi::GeometryFrozenMPI::GeometryFrozenMPI(int L_) {
    L = L_;
    V = L*L*L*L;
    neighbors_flat.resize(V*8);
    for (int t = 0; t < L; t++) {
        for (int z = 0; z < L; z++) {
            for (int y = 0; y < L; y++) {
                for (int x = 0; x < L; x++) {
                    size_t i = index(x, y, z, t);
                    neighbors_flat[index_neigh(i,0,0)] = index((x + 1) % L, y, z, t);
                    neighbors_flat[index_neigh(i,0,1)] = index((x - 1 + L) % L, y, z, t);
                    neighbors_flat[index_neigh(i,1,0)] = index(x, (y + 1) % L, z, t);
                    neighbors_flat[index_neigh(i,1,1)] = index(x, (y - 1 + L) % L, z, t);
                    neighbors_flat[index_neigh(i,2,0)] = index(x, y, (z + 1) % L, t);
                    neighbors_flat[index_neigh(i,2,1)] = index(x, y, (z - 1 + L) % L, t);
                    neighbors_flat[index_neigh(i,3,0)] = index(x, y, z, (t + 1) % L);
                    neighbors_flat[index_neigh(i,3,1)] = index(x, y, z, (t - 1 + L) % L);
                }
            }
        }
    }

    links_staples_flat.resize(V*3*6*4);
    for (int t = 0; t < L; t++) {
        for (int z = 0; z < L; z++) {
            for (int y = 0; y < L; y++) {
                for (int x = 0; x < L; x++) {
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
