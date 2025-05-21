#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aashen@berkeley.edu
#SBATCH --output=slurm_output/%x_%j.o
#SBATCH --error=slurm_output/%x_%j.e
#SBATCH -c 3
#SBATCH -p jsteinhardt

#––– Usage check
if [ -z "$1" ]; then
  echo "Usage: sbatch $0 {low|med|high}"
  exit 1
fi

#––– Validate argument
case "$1" in
  low|med|high) ;;
  *)
    echo "ERROR: Invalid argument '$1'"
    echo "  Must be one of: low, med, high"
    exit 1
    ;;
esac

#––– Dispatch
echo "→ Running sim1.R with level = $1"
Rscript sim1-divided.R "$1"
