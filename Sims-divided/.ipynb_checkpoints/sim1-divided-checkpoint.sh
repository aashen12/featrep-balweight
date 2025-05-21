#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aashen@berkeley.edu
#SBATCH --output=slurm_output/%x_%j.o
#SBATCH --error=slurm_output/%x_%j.e
#SBATCH -c 5
#SBATCH -p high

# —–––––––––––––––––––––––––––––––––––––––––––
# simple argument check
if [ -z "$1" ]; then
  echo "ERROR: you must pass one of {low, med, high}"
  echo "Usage: sbatch $0 {low|med|high}"
  exit 1
fi

case "$1" in
  low)   RSCRIPT="sim1-low.R"   ;;
  med)   RSCRIPT="sim1-med.R"   ;;
  high)  RSCRIPT="sim1-high.R"  ;;
  *)
    echo "Invalid argument: $1"
    echo "Usage: sbatch $0 {low|med|high}"
    exit 1
    ;;
esac

echo "→ Running Rscript $RSCRIPT"
Rscript "$RSCRIPT"
