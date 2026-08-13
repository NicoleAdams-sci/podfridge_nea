#!/bin/bash
#SBATCH --job-name=slim_pool
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:30:00
#SBATCH --output=logs/slim_pool_%j.out
#SBATCH --error=logs/slim_pool_%j.err

# ==============================================================================
# Slim an Unrelated Pool for the Focal Ranking Test — Shell Wrapper
# ==============================================================================
#
# Reads a pool CSV once and writes a small, fast-loading .fst copy filtered to
# the requested loci set(s) and the columns ranking needs. Run this ONCE per
# (pool, loci set); it is independent of pool generation, so changing the loci
# set only means re-running this, not regenerating the pool.
#
# Usage:
#   sbatch code/slim_unrelated_pool.sh <POOL_FILE> <LOCI_SETS> [OUTPUT_DIR]
#
# Example:
#   sbatch code/slim_unrelated_pool.sh \
#     output/unrelated_pool/all_N1000000_combined_unrelated_20260813_132702.csv \
#     core_13,expanded_20
#
# Local test:
#   bash code/slim_unrelated_pool.sh \
#     output/unrelated_pool/all_N1000_combined_unrelated_<dt>.csv core_13,expanded_20
# ==============================================================================

mkdir -p logs
module load Rtidyverse

if [ "$#" -lt 2 ]; then
    echo "Usage: sbatch code/slim_unrelated_pool.sh <POOL_FILE> <LOCI_SETS> [OUTPUT_DIR]" >&2
    echo "  e.g. sbatch code/slim_unrelated_pool.sh \\" >&2
    echo "         output/unrelated_pool/all_N1000000_combined_unrelated_<dt>.csv \\" >&2
    echo "         core_13,expanded_20" >&2
    exit 1
fi

echo "============================================================="
echo "  Slim Unrelated Pool"
echo "============================================================="
echo "  Arguments  : $@"
echo "  Started at : $(date)"
echo "============================================================="

Rscript code/slim_unrelated_pool.R "$@"

if [ $? -eq 0 ]; then
    echo "SUCCESS: Slim pool written at $(date)"
else
    echo "ERROR: Slim pool generation failed at $(date)"
    exit 1
fi
