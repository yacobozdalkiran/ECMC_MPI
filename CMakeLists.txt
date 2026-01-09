cmake_minimum_required(VERSION 4.1)
project(4d_MPI)

set(CMAKE_CXX_STANDARD 20)


find_package(MPI REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)

add_executable(4d_MPI main.cpp
        src/geometry/Geometry.cpp
        src/geometry/Geometry.h
        src/gauge/GaugeField.cpp
        src/gauge/GaugeField.h
        src/su3/utils.cpp
        src/su3/utils.h
        src/ecmc/ecmc.cpp
        src/ecmc/ecmc.h
        src/observables/observables.cpp
        src/observables/observables.h
        src/geometry/GeometryFrozen.cpp
        src/geometry/GeometryFrozen.h
        src/ecmc/ecmc_frozen.cpp
        src/ecmc/ecmc_frozen.h
        src/mpi/MpiTopology.cpp
        src/mpi/MpiTopology.h
        src/gauge/LocalGaugeField.cpp
        src/gauge/LocalGaugeField.h
        src/mpi/HalosExchange.cpp
        src/mpi/HalosExchange.h)

target_compile_options(4d_MPI PRIVATE
        -O3
        -ffast-math
        -march=native
        -Wall -Wextra
        -Wpedantic
        )

target_link_libraries(4d_MPI PRIVATE MPI::MPI_CXX Eigen3::Eigen)