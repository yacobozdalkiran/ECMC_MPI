cmake_minimum_required(VERSION 3.16)
project(ECMC_MPI)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(MPI REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)

if(NOT Eigen3_FOUND)
    message(STATUS "Eigen3 not found. Download via FetchContent...")

    include(FetchContent)
    FetchContent_Declare(
            eigen
            GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
            GIT_TAG 3.4.0
    )
    FetchContent_MakeAvailable(eigen)
else()
    message(STATUS "Eigen3 found on the system.")
endif()

add_executable(gauge_ecmc_mpi src/bin/gauge_ecmc_mpi.cpp
        src/gauge/GaugeField.cpp
        src/su3/utils.cpp
        src/observables/observables_mpi.cpp
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
            -O2
            #-xHost            # Remplace -march=native pour Intel
            -fp-model precise
            -Wall -Wextra
    )
else()
    target_compile_options(gauge_ecmc_mpi PRIVATE
            -O2
            #-ffast-math
            #-march=native
            -Wall -Wextra -Wpedantic
    )
endif()

target_link_libraries(gauge_ecmc_mpi PRIVATE MPI::MPI_CXX Eigen3::Eigen)

add_executable(gauge_single_core src/bin/gauge_single_core.cpp
        src/gauge/GaugeField.cpp
        src/geometry/Geometry.cpp
        src/ecmc/ecmc.cpp
        src/su3/utils.cpp
        src/observables/observables.cpp
        src/io/io.cpp
        src/io/params.h
)

target_link_libraries(gauge_single_core PRIVATE Eigen3::Eigen)

if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    target_compile_options(gauge_single_core PRIVATE
            -O2
            #-xHost
            -fp-model precise
            -Wall -Wextra
    )
else()
    target_compile_options(gauge_single_core PRIVATE
            -O2
            #-ffast-math
            #-march=native
            -Wall -Wextra -Wpedantic
    )
endif()

add_executable(gauge_metropolis src/bin/gauge_metropolis.cpp
        src/metropolis/metropolis.cpp
        src/gauge/GaugeField.cpp
        src/geometry/Geometry.cpp
        src/su3/utils.cpp
        src/observables/observables.cpp
        src/io/io.cpp
        src/io/params.h
)

target_link_libraries(gauge_metropolis PRIVATE Eigen3::Eigen)

if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    target_compile_options(gauge_metropolis PRIVATE
            -O2
            #-xHost
            -fp-model precise
            -Wall -Wextra
    )
else()
    target_compile_options(gauge_metropolis PRIVATE
            -O2
            #-ffast-math
            #-march=native
            -Wall -Wextra -Wpedantic
    )
endif()

add_executable(gauge_heatbath
        src/bin/gauge_heatbath.cpp
        src/heatbath/heatbath.cpp
        src/gauge/GaugeField.cpp
        src/geometry/Geometry.cpp
        src/su3/utils.cpp
        src/observables/observables.cpp
        src/io/io.cpp
        src/io/params.h
)

target_link_libraries(gauge_heatbath PRIVATE Eigen3::Eigen)

if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    target_compile_options(gauge_heatbath PRIVATE
            -O2
            #-xHost
            -fp-model precise
            -Wall -Wextra
)
else()
    target_compile_options(gauge_heatbath PRIVATE
            -O2
            #-ffast-math
            #-march=native
            -Wall -Wextra -Wpedantic
)
endif()
