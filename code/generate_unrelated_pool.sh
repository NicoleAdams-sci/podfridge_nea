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
# Generate Unrelated Pool / Database — Shell Wrapper
# ==============================================================================
#
# Usage:
#   sbatch code/generate_unrelated_pool.sh all 100000
#   sbatch code/generate_unrelated_pool.sh single 100000 Cauc
#   sbatch code/generate_unrelated_pool.sh mixed-proportions 100000 AfAm=0.15,Cauc=0.55,Hispanic=0.20,Asian=0.10
#   sbatch code/generate_unrelated_pool.sh mixed-counts AfAm=15000,Cauc=55000,Hispanic=20000,Asian=10000
#
# Local test:
#   bash code/generate_unrelated_pool.sh all 1000
#   bash code/generate_unrelated_pool.sh single 1000 Cauc
#   bash code/generate_unrelated_pool.sh mixed-proportions 1000 AfAm=0.15,Cauc=0.55,Hispanic=0.20,Asian=0.10
# ==============================================================================

mkdir -p logs
mkdir -p output/unrelated_pool

module load Rtidyverse

# Default behavior if no arguments supplied.
# You can change this default if desired.
if [ "$#" -eq 0 ]; then
    set -- all 250
fi

echo "============================================================="
echo "  Generate Unrelated Database"
echo "============================================================="
echo "  Arguments        : $@"
echo "  Output directory : output/unrelated_pool/"
echo "  Started at       : $(date)"
echo "============================================================="

Rscript code/generate_unrelated_pool.R "$@"

if [ $? -eq 0 ]; then
    echo "SUCCESS: Unrelated database generated at $(date)"
else
    echo "ERROR: Unrelated database generation failed at $(date)"
    exit 1
fi