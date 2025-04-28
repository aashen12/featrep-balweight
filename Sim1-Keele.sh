#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aashen@berkeley.edu
#SBATCH --error=slurm_output/%x_%j.e   # Error file
#SBATCH -c 5
#SBATCH -p yss

# Run the R script
Rscript Sim1-Keele.R