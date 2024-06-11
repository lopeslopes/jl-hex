#!/bin/bash

#SBATCH --nodes 1
#SBATCH --job-name=tbg_julia
#SBATCH --time=2-12
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
srun julia large_lattice.jl --oversubscribe
