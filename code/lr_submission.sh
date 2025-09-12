#!/bin/bash
# LR Calculation Job Submission Script
# This script generates the file list and submits the SLURM array job

echo "=== LR Calculation Job Submission ==="
echo "Started at: $(date)"

# Step 1: Generate file list of all pairs CSV files
echo "Generating file list..."
ls output/pairs_*_chunk*.csv > output/lr_file_list.txt 2>/dev/null

# Check if any files were found
if [ ! -s output/lr_file_list.txt ]; then
    echo "ERROR: No pairs_*.csv files found in output/ directory"
    echo "Make sure you have run the pair simulation jobs first"
    exit 1
fi

# Count files
FILE_COUNT=$(wc -l < output/lr_file_list.txt)
echo "Found $FILE_COUNT pairs files to process:"

# Show first few files as preview
head -5 output/lr_file_list.txt
if [ $FILE_COUNT -gt 5 ]; then
    echo "... and $(($FILE_COUNT - 5)) more files"
fi

# Step 2: Update the SLURM script with correct array count
echo "Updating SLURM script for $FILE_COUNT array tasks..."
sed "s/FILE_COUNT_PLACEHOLDER/$FILE_COUNT/" code/lr_wrapper.sh > code/lr_wrapper_ready.sh

# Step 3: Submit the job
echo "Submitting SLURM array job..."
JOB_ID=$(sbatch code/lr_wrapper_ready.sh | cut -d' ' -f4)

if [ $? -eq 0 ]; then
    echo "SUCCESS: Job submitted with ID: $JOB_ID"
    echo "Monitor with: squeue -u \$USER"
    echo "Check efficiency with: seff $JOB_ID"
    echo "View logs in: logs/lr_calc_${JOB_ID}_*.out"
    echo ""
    echo "LR results will be written to: output/LR/"
else
    echo "ERROR: Job submission failed"
    exit 1
fi

echo "=== Submission Complete ==="