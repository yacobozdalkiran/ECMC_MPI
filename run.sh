#!/bin/bash

#SBATCH --job-name=hot30000
#SBATCH --output=%x.o
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --partition=cpu_med

# Run MPI script
srun build/gauge_ecmc_mpi inputs/hot3000.txt