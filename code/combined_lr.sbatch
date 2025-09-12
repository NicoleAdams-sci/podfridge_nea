#!/bin/bash
#SBATCH --job-name=combined_lr
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=1G
#SBATCH --time=01:00:00
#SBATCH --array=1-FILE_COUNT_PLACEHOLDER
#SBATCH --output=logs/combined_lr_%A_%a.out
#SBATCH --error=logs/combined_lr_%A_%a.err

# Create logs directory if it doesn't exist
mkdir -p logs

# Create Combined_LR output directory if it doesn't exist
mkdir -p output/combined_LR

# Load required modules
module load Rtidyverse

# Read the specific LR file for this array task
LR_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" output/combined_lr_file_list.txt)

# Validate that we got a file
if [ -z "$LR_FILE" ]; then
    echo "ERROR: No file found for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

# Check if file exists
if [ ! -f "$LR_FILE" ]; then
    echo "ERROR: LR file does not exist: $LR_FILE"
    exit 1
fi

echo "Processing LR file: $LR_FILE"
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Started at: $(date)"

# Run the Combined LR calculation wrapper
Rscript code/combined_lr_wrapper.R "$LR_FILE"

# Check if the R script succeeded
if [ $? -eq 0 ]; then
    echo "SUCCESS: Combined LR calculation completed for $LR_FILE"
else
    echo "ERROR: Combined LR calculation failed for $LR_FILE"
    exit 1
fi

echo "Completed at: $(date)"
