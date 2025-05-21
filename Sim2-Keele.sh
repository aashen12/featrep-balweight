#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aashen@berkeley.edu
#SBATCH --output=slurm_output/%x_%j.o  # Output file
#SBATCH --error=slurm_output/%x_%j.e   # Error file
#SBATCH -c 10
#SBATCH -p yugroup

# Run the R script
Rscript Sim2-Keele.R