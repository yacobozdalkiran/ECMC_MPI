cmake_minimum_required(VERSION 3.16)
project(ECMC_MPI)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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
        src/io/io.cpp
        src/io/params.h
)


# Configuration des options selon le compilateur détecté
if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    # Options spécifiques à icpx (OneAPI)
    target_compile_options(gauge_ecmc_mpi PRIVATE
            -O3
            -xHost            # Remplace -march=native pour Intel
            -fp-model fast    # Équivalent plus robuste de -ffast-math
            -Wall -Wextra
    )
else()
    # On garde vos options pour GCC au cas où vous changeriez d'avis
    target_compile_options(gauge_ecmc_mpi PRIVATE
            -O3
            -ffast-math
            -march=native
            -Wall -Wextra -Wpedantic
    )
endif()

target_link_libraries(gauge_ecmc_mpi PRIVATE MPI::MPI_CXX Eigen3::Eigen)
