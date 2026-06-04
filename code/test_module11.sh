#!/bin/bash
#SBATCH --job-name=test_mod11
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --array=1-10
#SBATCH --output=logs/test_mod11_%A_%a.out
#SBATCH --error=logs/test_mod11_%A_%a.err

# ==============================================================================
# Module 11 Test — SLURM Array Wrapper
# ==============================================================================
#
# Runs one focal individual (replicate) per array task.
# Results saved to output/focal_ranking_test/lr_focal_<id>_<relative_type>.csv
# After all tasks complete, combine with:
#   Rscript code/combine_ranking_lrs.R <RELATIVE_TYPE>
#
# Usage:
#   sbatch code/test_module11.sh <FOCAL_FAMILY_FILE> <RELATIVE_TYPE> [FOCAL_POPULATION] [TESTED_POPULATIONS]
#
# Examples:
#   sbatch code/test_module11.sh output/focal_ranking_test/Asian_focal_parent_child1_full_siblings1_20260604_123505.csv parent_child
#   sbatch code/test_module11.sh output/focal_ranking_test/Asian_focal_second_cousins1_20260604_170702.csv second_cousins
#   sbatch code/test_module11.sh output/focal_ranking_test/Asian_focal_second_cousins1_20260604_170702.csv second_cousins Asian all
#
# Arguments:
#   $1 : FOCAL_FAMILY_FILE  — path to focal family CSV (required)
#   $2 : RELATIVE_TYPE      — true relative type to test (required)
#   $3 : FOCAL_POPULATION   — population label (default: Asian)
#   $4 : TESTED_POPULATIONS — allele freq population for LR calc (default: all)
#
# Prerequisites:
#   output/unrelated_pool/ : run generate_unrelated_pool.sh first
# ==============================================================================

mkdir -p logs
mkdir -p output/focal_ranking_test

module load Rtidyverse

FOCAL_FAMILY_FILE=${1:?"ERROR: FOCAL_FAMILY_FILE is required. Usage: sbatch test_module11.sh <FOCAL_FAMILY_FILE> <RELATIVE_TYPE>"}
RELATIVE_TYPE=${2:?"ERROR: RELATIVE_TYPE is required. Usage: sbatch test_module11.sh <FOCAL_FAMILY_FILE> <RELATIVE_TYPE>"}
FOCAL_POPULATION=${3:-Asian}
TESTED_POPULATIONS=${4:-all}
FOCAL_ID=${SLURM_ARRAY_TASK_ID}

# Validate file exists
if [ ! -f "${FOCAL_FAMILY_FILE}" ]; then
    echo "ERROR: Focal family file not found: ${FOCAL_FAMILY_FILE}"
    exit 1
fi

echo "============================================================="
echo "  Module 11 Test — Replicate ${FOCAL_ID}"
echo "============================================================="
echo "  Focal family file  : ${FOCAL_FAMILY_FILE}"
echo "  Relative type      : ${RELATIVE_TYPE}"
echo "  Focal population   : ${FOCAL_POPULATION}"
echo "  Tested populations : ${TESTED_POPULATIONS}"
echo "  Focal ID           : ${FOCAL_ID}"
echo "  Started at         : $(date)"
echo "============================================================="

Rscript code/test_module11.R \
    ${FOCAL_ID} \
    "${FOCAL_FAMILY_FILE}" \
    ${FOCAL_POPULATION} \
    ${TESTED_POPULATIONS} \
    ${RELATIVE_TYPE}

if [ $? -eq 0 ]; then
    echo "SUCCESS: Replicate ${FOCAL_ID} completed at $(date)"
else
    echo "ERROR: Replicate ${FOCAL_ID} failed at $(date)"
    exit 1
fi
