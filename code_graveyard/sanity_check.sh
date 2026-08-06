#!/bin/bash
#SBATCH --job-name=sanity_check
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=96G
#SBATCH --time=02:00:00
#SBATCH --output=logs/sanity_check_%j.out
#SBATCH --error=logs/sanity_check_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=nicolead@umich.edu

# SLURM wrapper for sanity_check_analysis.R + sanity_check_figures.R
# Usage: sbatch code/sanity_check.sh [OUTPUT_DIR] [COMBINED_LR_DIR]
# Example: sbatch code/sanity_check.sh
#          sbatch code/sanity_check.sh output/sanity_check output/combined_LR

mkdir -p logs

module load Rtidyverse

echo "================================================"
echo "Sanity Check Analysis"
echo "================================================"

OUTPUT_DIR=${1:-"output/sanity_check"}
COMBINED_LR_DIR=${2:-"output/combined_LR"}

COMBINED_COUNT=$(ls ${COMBINED_LR_DIR}/combined_LR_*.csv 2>/dev/null | wc -l)
echo "Found ${COMBINED_COUNT} combined_LR files in ${COMBINED_LR_DIR}"
echo ""

# Step 1: Run analysis
echo "Step 1: Running sanity check analysis..."
Rscript code/sanity_check_analysis.R "${OUTPUT_DIR}" "${COMBINED_LR_DIR}"
EXITCODE=$?

if [ ${EXITCODE} -ne 0 ]; then
    echo "ERROR: Analysis failed with exit code ${EXITCODE}"
    exit ${EXITCODE}
fi

# Step 2: Generate figures
echo ""
echo "Step 2: Generating figures..."
Rscript code/sanity_check_figures.R "${OUTPUT_DIR}"
EXITCODE=$?

if [ ${EXITCODE} -ne 0 ]; then
    echo "WARNING: Figure generation failed (exit code ${EXITCODE})"
    echo "Analysis CSVs were still produced successfully."
fi

echo ""
echo "================================================"
echo "Output files in: ${OUTPUT_DIR}"
ls -la "${OUTPUT_DIR}/"
echo ""
echo "Job completed at: $(date)"
echo "Total runtime: ${SECONDS} seconds"
echo "================================================"

seff ${SLURM_JOB_ID}
