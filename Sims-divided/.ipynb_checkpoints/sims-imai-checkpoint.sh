#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aashen@berkeley.edu
#SBATCH --output=slurm_output/%x_%j.o
#SBATCH --error=slurm_output/%x_%j.e
#SBATCH -c 22
#SBATCH -p jsteinhardt

Rscript sims-imai.R
