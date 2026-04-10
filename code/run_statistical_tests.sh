#!/bin/bash
#SBATCH --job-name=stats_tests
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/stats_tests_%j.out
#SBATCH --error=logs/stats_tests_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=nicolead@umich.edu

# SLURM wrapper for run_statistical_tests.R
# Loads combined_LR_all.rds for Section 2 (mismatched population pair-level tests).
# Sections 1 and 3 read smaller files but are included in the same job for
# a single consolidated run.
#
# Usage:   sbatch code/run_statistical_tests.sh <input_dir>
# Example: sbatch code/run_statistical_tests.sh output/lr_analysis_20260410

mkdir -p logs

module load Rtidyverse

if [ -z "$1" ]; then
    echo "ERROR: Input directory required"
    echo "Usage: sbatch code/run_statistical_tests.sh <input_dir>"
    exit 1
fi

INPUT_DIR=$1

if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory does not exist: $INPUT_DIR"
    exit 1
fi

# Check all three required input files exist before starting
MISSING=0
for f in combined_LR_all.rds combined_LR_match.csv.gz proportions_with_classification.csv; do
    if [ ! -f "$INPUT_DIR/$f" ]; then
        echo "ERROR: Required file not found: $INPUT_DIR/$f"
        MISSING=1
    fi
done

if [ $MISSING -eq 1 ]; then
    echo ""
    echo "Run the following first:"
    echo "  sbatch code/analyze_lr_outputs.sh $INPUT_DIR"
    echo "  sbatch code/prepare_combined_lr_intermediates.sh $INPUT_DIR"
    exit 1
fi

echo "================================================"
echo "Statistical Tests Pipeline"
echo "================================================"
echo "Input directory: $INPUT_DIR"
echo ""

Rscript code/run_statistical_tests.R "$INPUT_DIR"
EXITCODE=$?

echo "================================================"

if [ ${EXITCODE} -eq 0 ]; then
    echo "SUCCESS: Statistical tests complete"
    echo ""
    echo "Output files:"
    for f in stats_matched.csv stats_mismatched_population.csv stats_fpr_cutoffs.csv; do
        if [ -f "$INPUT_DIR/stats/$f" ]; then
            size=$(du -h "$INPUT_DIR/stats/$f" | cut -f1)
            echo "  - stats/$f ($size)"
        fi
    done
else
    echo "ERROR: run_statistical_tests.R failed with exit code ${EXITCODE}"
    echo "Check logs/stats_tests_${SLURM_JOB_ID}.out for details"
    exit ${EXITCODE}
fi

echo ""
echo "================================================"
echo "Job completed at: $(date)"
echo "Total runtime: ${SECONDS} seconds"
echo "================================================"
