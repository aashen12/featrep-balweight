#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aashen@berkeley.edu
#SBATCH --output=slurm_output/%x_%j.o
#SBATCH --error=slurm_output/%x_%j.e
#SBATCH -c 3
#SBATCH -p jsteinhardt

# allowed values
PARAMS=(1 2.5 5 7.5 10)

# 1) make sure we got an argument
if [ -z "$1" ]; then
  echo "Usage: sbatch $0 <overlap_degree>"
  echo "  must be one of: ${PARAMS[*]}"
  exit 1
fi

# 2) check it’s in our list
VALID=false
for v in "${PARAMS[@]}"; do
  if [[ "$1" == "$v" ]]; then
    VALID=true
    break
  fi
done

if ! $VALID; then
  echo "ERROR: invalid overlap_degree '$1'"
  echo "  must be one of: ${PARAMS[*]}"
  exit 1
fi

echo "→ Launching sim1.R with overlap_degree = $1"
Rscript sim1-single.R "$1"
