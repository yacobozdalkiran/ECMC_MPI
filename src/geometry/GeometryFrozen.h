//
// Created by ozdalkiran-l on 1/8/26.
//

#ifndef INC_4D_MPI_GEOMETRYFROZEN_H
#define INC_4D_MPI_GEOMETRYFROZEN_H

#include <vector>

//Geometry of a lattice with periodic BC
class GeometryFrozen {
public:
    int L;
    int T;
    size_t V;
private:
    std::vector<size_t> neighbors_flat;
    std::vector<std::pair<size_t, int>> links_staples_flat;
    std::vector<bool> is_frozen_flat;

public:
    [[nodiscard]] size_t index(int x, int y, int z, int t) const {
        return ((static_cast<size_t>(t)*L + z) * L + y)*L + x;
    }

    [[nodiscard]] static size_t index_neigh(size_t site, int mu, int dir) {
        return site*8 + mu*2 + dir;
    }

    [[nodiscard]] static size_t index_staples(size_t site, int mu, int i_staple, int i_link) {
        return site*3*6*4 + mu*3*6 + i_staple*3 + i_link;
    }

    [[nodiscard]] static size_t index_frozen(size_t site, int mu) {
        return site*4 + mu;
    }

    explicit GeometryFrozen(int L, int T);

    //Returns the neighbor of site in direction +mu if dir=0, -mu if dir=1
    [[nodiscard]] size_t get_neigh(size_t site, int mu, int dir) const {
        return neighbors_flat[index_neigh(site, mu, dir)];
    }

    //Returns the coords of the link of index i_link of the staple of index i_staple of <site,mu>
    [[nodiscard]] std::pair<size_t,int> get_link_staple(size_t site, int mu, int i_staple, int i_link) const {
        return links_staples_flat[index_staples(site, mu, i_staple, i_link)];
    }

    [[nodiscard]] bool is_frozen(size_t site, int mu) const {
        return is_frozen_flat[index_frozen(site, mu)];
    }
};

#endif //INC_4D_MPI_GEOMETRYFROZEN_H
