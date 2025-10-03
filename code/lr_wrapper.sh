#!/bin/bash
#SBATCH --job-name=lr_calc
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH --time=04:00:00
#SBATCH --array=1-FILE_COUNT_PLACEHOLDER
#SBATCH --output=logs/lr_calc_%A_%a.out
#SBATCH --error=logs/lr_calc_%A_%a.err

# Create logs directory if it doesn't exist
mkdir -p logs

# Create LR output directory if it doesn't exist
mkdir -p output/LR

# Load required modules
module load Rtidyverse

# Read the specific pairs file for this array task
PAIRS_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" output/lr_file_list.txt)

# Validate that we got a file
if [ -z "$PAIRS_FILE" ]; then
    echo "ERROR: No file found for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

# Check if file exists
if [ ! -f "$PAIRS_FILE" ]; then
    echo "ERROR: Pairs file does not exist: $PAIRS_FILE"
    exit 1
fi

echo "Processing pairs file: $PAIRS_FILE"
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Started at: $(date)"

# Run the LR calculation wrapper
Rscript code/lr_wrapper.R "$PAIRS_FILE"

# Check if the R script succeeded
if [ $? -eq 0 ]; then
    echo "SUCCESS: LR calculation completed for $PAIRS_FILE"
else
    echo "ERROR: LR calculation failed for $PAIRS_FILE"
    exit 1
fi

echo "Completed at: $(date)"