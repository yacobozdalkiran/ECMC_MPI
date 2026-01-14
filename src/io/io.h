//
// Created by ozdalkiran-l on 1/14/26.
//

#ifndef ECMC_MPI_IO_H
#define ECMC_MPI_IO_H

#include <vector>

#include "params.h"

namespace io {
    //Output
    void save_double(const std::vector<double> &data, const std::string &filename, int precision);
    void save_double_params(const std::vector<double> &data, const RunParams &params, const std::string &filename, int precision);

    //Input
    std::string trim(const std::string& s);
    void load_params(const std::string& filename, RunParams& rp);
}


#endif //ECMC_MPI_IO_H