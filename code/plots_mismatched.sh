#!/bin/bash
#SBATCH --job-name=plot_mismatch
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=21G
#SBATCH --time=01:00:00
#SBATCH --output=logs/plot_mismatch_pop_%j.out
#SBATCH --error=logs/plot_mismatch_pop_%j.err

# Usage: sbatch code/plot_mismatched_wrapper.sh <input_dir> [output_subdir]
# Example: sbatch code/plot_mismatched_wrapper.sh output/lr_analysis_20250105
# Example with custom output: sbatch code/plot_mismatched_wrapper.sh output/lr_analysis_20250105 mismatched_plots_final

# Create logs directory if it doesn't exist
mkdir -p logs

# Load required modules
module load Rtidyverse

# Get input directory from command line argument or use default
if [ -z "$1" ]; then
    echo "ERROR: Input directory required"
    echo "Usage: sbatch code/plot_mismatched_wrapper.sh <input_dir> [output_subdir]"
    exit 1
fi

INPUT_DIR=$1

# Check if input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory does not exist: $INPUT_DIR"
    exit 1
fi

# Get optional output subdirectory name
OUTPUT_SUBDIR=${2:-""}


echo "Input directory: $INPUT_DIR"

# Run the plotting script
if [ -n "$OUTPUT_SUBDIR" ]; then
    Rscript code/plots_mismatched.R "$INPUT_DIR" "$OUTPUT_SUBDIR"
else
    Rscript code/plots_mismatched.R "$INPUT_DIR"
fi