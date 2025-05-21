#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aashen@berkeley.edu
#SBATCH --output=slurm_output/%x_%j.o
#SBATCH --error=slurm_output/%x_%j.e
#SBATCH -c 3
#SBATCH -p high

# check for an argument
if [ -z "$1" ]; then
  echo "Usage: sbatch $0 <param>"
  echo "  <param> must be one of: 1, 2.5, 5, 7.5"
  exit 1
fi

# remove the decimal point to match your filenames:
#   2.5 → 25, 7.5 → 75, 1 → 1, 5 → 5
FILE_SUFFIX="${1//./}"

RSCRIPT="sim1-${FILE_SUFFIX}.R"

# sanity check
if [ ! -f "$RSCRIPT" ]; then
  echo "ERROR: $RSCRIPT not found!"
  exit 1
fi

echo "→ Running Rscript $RSCRIPT"
Rscript "$RSCRIPT"
