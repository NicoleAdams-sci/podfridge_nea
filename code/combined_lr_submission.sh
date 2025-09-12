#!/bin/bash
# Combined LR Calculation Job Submission Script
# Usage: code/combined_lr_submission.sh [file_list.txt]

# Use provided file list or create default in output directory
FILE_LIST=${1:-output/combined_lr_file_list.txt}

# If using default, generate it from LR files
if [ "$FILE_LIST" = "output/combined_lr_file_list.txt" ]; then
    ls output/LR/LR_*.csv > "$FILE_LIST"
fi

# Count files and update SLURM script with file count AND file list path
FILE_COUNT=$(wc -l < "$FILE_LIST")
sed -e "s/FILE_COUNT_PLACEHOLDER/$FILE_COUNT/" \
    -e "s|FILE_LIST_PLACEHOLDER|$FILE_LIST|" \
    code/combined_lr.sh > code/combined_lr_ready.sh

# Submit job
sbatch code/combined_lr_ready.sh