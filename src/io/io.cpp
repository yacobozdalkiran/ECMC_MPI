//
// Created by ozdalkiran-l on 1/14/26.
//

#include <fstream>
#include <iostream>
#include <filesystem>
#include <mpi.h>
#include <map>
#include <string>
#include <sstream>


#include "io.h"

namespace fs = std::filesystem;

//Saves a vector of doubles in ../data/filename.txt
void io::save_double(const std::vector<double> &data, const std::string &filename, int precision) {
    //Create a data folder if doesn't exists
    fs::path dir("../data");

    try {
        if (!fs::exists(dir)) {
            fs::create_directories(dir);
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Couldn't create folder data : " << e.what() << std::endl;
        return;
    }
    fs::path filepath = dir/ (filename+".txt");

    int counter = 1;
    while (fs::exists(filepath)) {
        //new name generation if file already exists
        std::string new_name = filename + "_" + std::to_string(counter) + ".txt";
        filepath = dir / new_name;
        counter++;
    }



    std::ofstream file(filepath);
    if (!file.is_open()) {
        std::cout << "Could not open file " << filepath << "\n";
        return;
    }
    file << std::setprecision(precision);
    for (const double &x : data) {
        file << x << "\n";
    }
    file.close();
    std::cout << "Data written in " << filepath << "\n";
}

//Saves a vector of doubles in ../data/filename and the run parameters in ../data/filename_params
void io::save_double_params(const std::vector<double> &data, const RunParams &params, const std::string &filename,
    int precision) {
    io::save_double(data, filename, precision);
    //No need to create data folder
    fs::path dir("../data");

    //We create the params file
    std::string param_filename = filename + "_params.txt";
    fs::path filepath = dir/param_filename;

    int counter = 1;
    while (fs::exists(filepath)) {
        //new name generation if file already exists
        std::string new_name = filename + "_" + std::to_string(counter) + "_params.txt";
        filepath = dir / new_name;
        counter++;
    }

    std::ofstream file(filepath);
    if (!file.is_open()) {
        std::cout << "Could not open file " << filepath << "\n";
        return;
    }
    //We write the run params
    file << "L_tot: " << params.L_core*params.n_core_dims << "\n"
            << "L_core: " << params.L_core << "\n"
            << "n_core_dims: " << params.n_core_dims << "\n"
            << "L_shift: " << params.L_shift << "\n"
            << "n_shift: " << params.n_shift << "\n"
            << "beta: " << params.ecmc_params.beta << "\n"
            << "N_samples: " << params.ecmc_params.N_samples << "\n"
            << "param_theta_sample: " << params.ecmc_params.param_theta_sample << "\n"
            << "param_theta_refresh: " << params.ecmc_params.param_theta_refresh << "\n"
            << "poisson: " << (params.ecmc_params.poisson ? "true" : "false") << "\n"
            << "epsilon_set: " << params.ecmc_params.epsilon_set;
    file.close();
}

//Utilitary function to trim the spaces
std::string io::trim(const std::string& s) {
    size_t first = s.find_first_not_of(" \t");
    if (first == std::string::npos) return "";
    size_t last = s.find_last_not_of(" \t");
    return s.substr(first, (last - first + 1));
}

//Loads the parameters contained in filename into RunParams rp
void io::load_params(const std::string& filename, RunParams& rp) {
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Impossible d'ouvrir " + filename);

    std::map<std::string, std::string> config;
    std::string line;

    while (std::getline(file, line)) {
        //Ignore comments and empty lines
        if (line.empty() || line[0] == '#') continue;

        std::istringstream is_line(line);
        std::string key, value;
        if (std::getline(is_line, key, '=') && std::getline(is_line, value)) {
            config[trim(key)] = trim(value);
        }
    }

    //Lattice params
    if (config.count("L_core"))      rp.L_core = std::stoi(config["L_core"]);
    if (config.count("n_core_dims")) rp.n_core_dims = std::stoi(config["n_core_dims"]);
    if (config.count("cold_start"))  rp.cold_start = (config["cold_start"] == "true");
    if (config.count("L_shift"))     rp.L_shift = std::stoi(config["L_shift"]);
    if (config.count("n_shift"))     rp.n_shift = std::stoi(config["n_shift"]);
    if (config.count("stype_pos"))   rp.stype_pos = (config["stype_pos"] == "true");

    //ECMC params
    if (config.count("beta"))                rp.ecmc_params.beta = std::stod(config["beta"]);
    if (config.count("N_samples"))           rp.ecmc_params.N_samples = std::stoi(config["N_samples"]);
    if (config.count("param_theta_sample"))  rp.ecmc_params.param_theta_sample = std::stod(config["param_theta_sample"]);
    if (config.count("param_theta_refresh")) rp.ecmc_params.param_theta_refresh = std::stod(config["param_theta_refresh"]);
    if (config.count("poisson"))             rp.ecmc_params.poisson = (config["poisson"] == "true");
    if (config.count("epsilon_set"))         rp.ecmc_params.epsilon_set = std::stod(config["epsilon_set"]);
}