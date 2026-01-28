#!/bin/bash

#SBATCH --job-name=hbhot30k
#SBATCH --output=%x.o
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --partition=cpu_med

#Source necessary modules
source modules_load.sh

# Run MPI script
srun build/gauge_ecmc_mpi inputs/heatbath_mpi_benchmark/hot30k.txt