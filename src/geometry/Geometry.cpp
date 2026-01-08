//
// Created by ozdalkiran-l on 1/8/26.
//

#include "Geometry.h"

Geometry::Geometry(int L_, int T_) {
    L = L_;
    T = T_;
    V = T*L*L*L;
    neighbors_flat.resize(V*8);
    for (int t = 0; t < T; t++) {
        for (int z = 0; z < L; z++) {
            for (int y = 0; y < L; y++) {
                for (int x = 0; x < L; x++) {
                    size_t i = index(x, y, z, t);
                    neighbors_flat[i*8 + 0*2 + 0] = index((x + 1) % L, y, z, t);
                    neighbors_flat[i*8 + 0*2 + 1] = index((x - 1 + L) % L, y, z, t);
                    neighbors_flat[i*8 + 1*2 + 0] = index(x, (y + 1) % L, z, t);
                    neighbors_flat[i*8 + 1*2 + 1] = index(x, (y - 1 + L) % L, z, t);
                    neighbors_flat[i*8 + 2*2 + 0] = index(x, y, (z + 1) % L, t);
                    neighbors_flat[i*8 + 2*2 + 1] = index(x, y, (z - 1 + L) % L, t);
                    neighbors_flat[i*8 + 3*2 + 0] = index(x, y, z, (t + 1) % T);
                    neighbors_flat[i*8 + 3*2 + 1] = index(x, y, z, (t - 1 + T) % T);
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

                            links_staples_flat[site*3*6*4+ mu*3*6 + j*3 + 0] = {xmu, nu};
                            links_staples_flat[site*3*6*4+ mu*3*6 + j*3 + 1] = {xnu, mu};
                            links_staples_flat[site*3*6*4+ mu*3*6 + j*3 + 2] = {site, nu};

                            links_staples_flat[site*3*6*4+ mu*3*6 + (j+1)*3 + 0] = {xmunu, nu};
                            links_staples_flat[site*3*6*4+ mu*3*6 + (j+1)*3 + 1] = {xmnu, mu};
                            links_staples_flat[site*3*6*4+ mu*3*6 + (j+1)*3 + 2] = {xmnu, nu};

                            j += 2;
                        }
                    }
                }
            }
        }
    }
}
