#!/bin/bash
#SBATCH --job-name=unrelated_pool
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/unrelated_pool_%j.out
#SBATCH --error=logs/unrelated_pool_%j.err

# ==============================================================================
# Generate Unrelated Pool — Shell Wrapper
# ==============================================================================
#
# Generates a fixed pool of unrelated individuals across all four populations,
# saved to output/unrelated_pool/. Run once and reuse across all ranking tests.
#
# Usage:
#   sbatch code/generate_unrelated_pool.sh                     # all defaults
#   sbatch code/generate_unrelated_pool.sh 250                 # 250 per pop (test)
#   sbatch code/generate_unrelated_pool.sh 100000              # full run
#   bash code/generate_unrelated_pool.sh 250                   # local test
#
# Arguments (optional):
#   $1 : N_UNRELATED_PER_POP — unrelated individuals per population (default: 250)
#
# Output:
#   output/unrelated_pool/<population>_unrelated_<datetime>.csv (one per population)
# ==============================================================================

mkdir -p logs
mkdir -p output/unrelated_pool

module load Rtidyverse

N_UNRELATED_PER_POP=${1:-250}

echo "============================================================="
echo "  Generate Unrelated Pool"
echo "============================================================="
echo "  Unrelated per population : ${N_UNRELATED_PER_POP}"
echo "  Total individuals        : $(( N_UNRELATED_PER_POP * 4 ))"
echo "  Output directory         : output/unrelated_pool/"
echo "  Started at               : $(date)"
echo "============================================================="

Rscript code/generate_unrelated_pool.R ${N_UNRELATED_PER_POP}

if [ $? -eq 0 ]; then
    echo "SUCCESS: Unrelated pool generated at $(date)"
else
    echo "ERROR: Unrelated pool generation failed at $(date)"
    exit 1
fi
