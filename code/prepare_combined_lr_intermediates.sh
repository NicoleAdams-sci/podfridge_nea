#!/bin/bash
#SBATCH --job-name=prep_intermediates
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=48G
#SBATCH --time=01:00:00
#SBATCH --output=logs/prep_intermediates_%j.out
#SBATCH --error=logs/prep_intermediates_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=nicolead@umich.edu

# SLURM wrapper for prepare_combined_lr_intermediates.R
# Loads combined_LR_all.rds and writes three aggregated CSVs for plotting.
#
# Usage:   sbatch code/prepare_combined_lr_intermediates.sh <input_dir>
# Example: sbatch code/prepare_combined_lr_intermediates.sh output/lr_analysis_20260130

mkdir -p logs

module load Rtidyverse

if [ -z "$1" ]; then
    echo "ERROR: Input directory required"
    echo "Usage: sbatch code/prepare_combined_lr_intermediates.sh <input_dir>"
    exit 1
fi

INPUT_DIR=$1

if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory does not exist: $INPUT_DIR"
    exit 1
fi

if [ ! -f "$INPUT_DIR/combined_LR_all.rds" ]; then
    echo "ERROR: combined_LR_all.rds not found in $INPUT_DIR"
    echo "Run analyze_lr_outputs.sh first."
    exit 1
fi

echo "================================================"
echo "Prepare Intermediates Pipeline"
echo "================================================"
echo "Input directory: $INPUT_DIR"
echo ""

Rscript code/prepare_combined_lr_intermediates.R "$INPUT_DIR"
EXITCODE=$?

echo "================================================"

if [ ${EXITCODE} -eq 0 ]; then
    echo "SUCCESS: Intermediates written to $INPUT_DIR"
    echo ""
    echo "Output files:"
    for f in proportions_with_classification.csv mismatched_pop_robustness.csv mismatched_pop_heatmap.csv; do
        if [ -f "$INPUT_DIR/$f" ]; then
            size=$(du -h "$INPUT_DIR/$f" | cut -f1)
            echo "  - $f ($size)"
        fi
    done
    echo ""
    echo "Next steps — run plotting scripts interactively:"
    echo "  Rscript code/plots_cutoffs_publication.R $INPUT_DIR"
    echo "  Rscript code/plots_mismatched_population.R $INPUT_DIR"
    echo "  Rscript code/plots_matched_publication.R $INPUT_DIR"
else
    echo "ERROR: prepare_combined_lr_intermediates failed with exit code ${EXITCODE}"
    echo "Check logs/prep_intermediates_${SLURM_JOB_ID}.out for details"
    exit ${EXITCODE}
fi

echo ""
echo "================================================"
echo "Job completed at: $(date)"
echo "Total runtime: ${SECONDS} seconds"
echo "================================================"
