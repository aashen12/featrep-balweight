#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aashen@berkeley.edu
#SBATCH --output=slurm_output/%x_%j.o
#SBATCH --error=slurm_output/%x_%j.e
#SBATCH -c 26
#SBATCH -p jsteinhardt

PARAMS=(30 100)

overlap=$1

# 2) Validate overlap_degree
VALID=false
for v in "${PARAMS[@]}"; do
  if [[ "$overlap" == "$v" ]]; then
    VALID=true
    break
  fi
done

if ! $VALID; then
  echo "ERROR: invalid overlap_degree '$overlap'"
  echo "  must be one of: ${PARAMS[*]}"
  exit 1
fi



echo "→ Running sims-single.R with overlap_degree=$overlap"
Rscript sims-jose.R "$overlap"