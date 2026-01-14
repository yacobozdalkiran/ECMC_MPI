#!/bin/bash

#SBATCH --job-name=beta_array
#SBATCH --output=%x_%a_%A.o
#SBATCH --time=00:50:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --array=1-9
#SBATCH --partition=cpu_med

# Run MPI script
srun build/gauge_ecmc_mpi inputs/beta${SLURM_ARRAY_TASK_ID}.txt