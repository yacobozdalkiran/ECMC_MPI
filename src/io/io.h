//
// Created by ozdalkiran-l on 1/14/26.
//

#ifndef ECMC_MPI_IO_H
#define ECMC_MPI_IO_H

#include <vector>
#include <sstream>
#include <iomanip>

#include "params.h"

namespace io {
    //Output
    void save_double(const std::vector<double> &data, const std::string &filename, int precision);
    void save_double_params(const std::vector<double> &data, const RunParams &params, const std::string &filename, int precision);

    //Input
    std::string trim(const std::string& s);
    void load_params(const std::string& filename, RunParams& rp);
    void load_params_sc(const std::string &filename, RunParamsSC& rp);

    //Utilitaries
    inline std::string format_double(double val, int precision) {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(precision) << val;
        return ss.str();
    }
}


#endif //ECMC_MPI_IO_H