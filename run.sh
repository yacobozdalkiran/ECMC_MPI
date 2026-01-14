#!/bin/bash

#SBATCH --job-name=test_lat4d_mpi
#SBATCH --output=%x.o%j 
#SBATCH --time=00:20:00 
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --partition=cpu_short

# Run MPI script
srun build/gauge_ecmc_mpi inputs/test.txt