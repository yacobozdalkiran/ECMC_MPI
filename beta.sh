#!/bin/bash

#SBATCH --job-name=beta
#SBATCH --output=%x_%a_%A.o
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --array=1-9
#SBATCH --partition=cpu_med

#Source necessary modules
source modules_load.sh

# Run MPI script
srun build/gauge_ecmc_mpi inputs/beta${SLURM_ARRAY_TASK_ID}.txt