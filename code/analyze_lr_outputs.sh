#!/bin/bash
#SBATCH --job-name=analyze_lr
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=96G
#SBATCH --time=02:00:00
#SBATCH --output=logs/analyze_lr_%j.out
#SBATCH --error=logs/analyze_lr_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=nicolead@umich.edu

# SLURM wrapper for analyze_lr_outputs.R
# Usage: sbatch code/analyze_lr.sh [output_dir]
# Example: sbatch code/analyze_lr.sh
#          sbatch code/analyze_lr.sh output/lr_analysis_custom

# Create logs directory if it doesn't exist
mkdir -p logs

# Load required modules
module load Rtidyverse

echo "================================================"
echo "LR Analysis Pipeline"
echo "================================================"

# Get optional output directory from command line
OUTPUT_DIR=${1:-""}


# Count input files for verification
LR_COUNT=$(ls output/LR/LR_*.csv 2>/dev/null | wc -l)
COMBINED_COUNT=$(ls output/combined_LR/combined_LR_*.csv 2>/dev/null | wc -l)

echo "Found input files:"
echo "  - Single-locus LR files: ${LR_COUNT}"
echo "  - Combined LR files: ${COMBINED_COUNT}"
echo ""


# Run the analysis script
echo "Running LR analysis pipeline..."
echo "================================================"

if [ -z "${OUTPUT_DIR}" ]; then
    # Run without specifying output directory (will use default timestamp)
    Rscript code/analyze_lr_outputs.R
    EXITCODE=$?
else
    # Run with specified output directory
    Rscript code/analyze_lr_outputs.R "${OUTPUT_DIR}"
    EXITCODE=$?
fi

echo "================================================"

# Check if R script succeeded
if [ ${EXITCODE} -eq 0 ]; then
    echo "SUCCESS: LR analysis completed successfully"
    
    # Find the output directory that was created
    if [ -z "${OUTPUT_DIR}" ]; then
        # Find most recent lr_analysis directory
        OUTPUT_DIR=$(ls -dt output/lr_analysis_* 2>/dev/null | head -1)
    fi
    
    if [ -d "${OUTPUT_DIR}" ]; then
        echo "Output saved to: ${OUTPUT_DIR}"
        echo ""
        echo "Summary file counts:"
        for pop in AfAm Cauc Hispanic Asian; do
            if [ -d "${OUTPUT_DIR}/${pop}_summary" ]; then
                count=$(ls ${OUTPUT_DIR}/${pop}_summary/*.csv 2>/dev/null | wc -l)
                echo "  - ${pop}: ${count} summary files"
            fi
        done
        echo ""
        echo "Next steps:"
        echo "1. Run plotting scripts:"
        echo "   Rscript code/plots_known_NEA.R ${OUTPUT_DIR}"
        echo "   Rscript code/plots_known_vs_tested_NEA.R ${OUTPUT_DIR}"
        echo "   Rscript code/plots_known_vs_tested_byLocus_NEA.R ${OUTPUT_DIR}"
    fi
else
    echo "ERROR: LR analysis failed with exit code ${EXITCODE}"
    echo "Check the log files for details:"
    echo "  - Output log: logs/analyze_lr_${SLURM_JOB_ID}.out"
    echo "  - Error log: logs/analyze_lr_${SLURM_JOB_ID}.err"
    exit ${EXITCODE}
fi

echo ""
echo "================================================"
echo "Job completed at: $(date)"
echo "Total runtime: ${SECONDS} seconds"
echo "================================================"

# Generate efficiency report
seff ${SLURM_JOB_ID}