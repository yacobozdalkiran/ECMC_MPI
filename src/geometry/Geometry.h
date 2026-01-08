//
// Created by ozdalkiran-l on 1/8/26.
//

#ifndef INC_4D_MPI_GEOMETRY_H
#define INC_4D_MPI_GEOMETRY_H

#include <vector>
#include <array>

//Geometry of a lattice with periodic BC
struct Geometry {
    int L;
    int T;
    size_t V;
    std::vector<size_t> neighbors_flat;
    std::vector<std::pair<size_t, int>> links_staples_flat;

    [[nodiscard]] size_t index(int x, int y, int z, int t) const {
        return ((static_cast<size_t>(t)*L + z) * L + y)*L + x;
    }

    explicit Geometry(int L, int T);
    [[nodiscard]] size_t get_neigh(size_t site, int mu, int dir) const {
        return neighbors_flat[site*8 + mu*2 + dir];
    }

    //Returns the coords of the link of index i_link of the staple of index i_staple of <site,mu>
    [[nodiscard]] std::pair<size_t,int> get_link_staple(size_t site, int mu, int i_staple, int i_link) const {
        return links_staples_flat[site*3*6*4 + mu*3*6 + i_staple*3 + i_link];
    }
};


#endif //INC_4D_MPI_GEOMETRY_H