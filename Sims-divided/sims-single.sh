#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aashen@berkeley.edu
#SBATCH --output=slurm_output/%x_%j.o
#SBATCH --error=slurm_output/%x_%j.e
#SBATCH -c 5
#SBATCH -p high

# allowed numeric overlap degrees
PARAMS=(1 2.5 5 7.5 10)

# 1) Make sure we got exactly two args
if [ "$#" -ne 2 ]; then
  echo "Usage: sbatch $0 <overlap_degree> <sim_num>"
  echo "  <overlap_degree> ∈ ${PARAMS[*]}"
  echo "  <sim_num> must be a positive integer"
  exit 1
fi

overlap=$1
sim_num=$2

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

# 3) Validate sim_num is an integer
if ! [[ "$sim_num" =~ ^[0-9]+$ ]]; then
  echo "ERROR: sim_num must be a positive integer, got '$sim_num'"
  exit 1
fi

echo "→ Running sims-single.R with overlap_degree=$overlap, sim_num=$sim_num"
Rscript sims-single.R "$overlap" "$sim_num"
