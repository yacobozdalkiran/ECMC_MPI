cmake_minimum_required(VERSION 4.1)
project(ECMC_MPI)

set(CMAKE_CXX_STANDARD 20)


find_package(MPI REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)

add_executable(gauge_ecmc_mpi src/bin/gauge_ecmc_mpi.cpp
        src/gauge/GaugeField.cpp
        src/su3/utils.cpp
        src/observables/observables.cpp
        src/geometry/GeometryFrozen.cpp
        src/mpi/MpiTopology.cpp
        src/mpi/HalosExchange.cpp
        src/mpi/Halo.h
        src/ecmc/ecmc_mpi.cpp
        include/cxxopts.hpp
        src/io/io.cpp
        src/io/params.cpp
        src/io/params.h
)


target_compile_options(gauge_ecmc_mpi PRIVATE
        -O3
        -ffast-math
        -march=native
        -Wall -Wextra
        -Wpedantic
)

target_link_libraries(gauge_ecmc_mpi PRIVATE MPI::MPI_CXX Eigen3::Eigen)
