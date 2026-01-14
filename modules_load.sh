#!/bin/bash
module purge
module load gcc/13.2.0/gcc-4.8.5
# Charger les compilateurs OneAPI permet souvent d'exposer TBB
module load intel-oneapi-compilers/2023.2.1/gcc-11.2.0
module load intel-oneapi-mpi/2021.10.0/gcc-13.2.0
module load cmake/3.21.4/gcc-13.2.0