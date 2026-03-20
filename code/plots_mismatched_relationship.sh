#!/bin/bash
#SBATCH --job-name=plot_mismatch_relationship
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=42G
#SBATCH --time=01:00:00
#SBATCH --output=logs/plot_mismatch_relationship_%j.out
#SBATCH --error=logs/plot_mismatch_relationship_%j.err

# Usage: sbatch code/plots_mismatched_relationship.sh <input_dir> [output_subdir]
# Example: sbatch code/plots_mismatched_relationship.sh output/lr_analysis_20250105
# Example with custom output: sbatch code/plots_mismatched_relationship.sh output/lr_analysis_20250105 output/lr_analysis_20250105/plots_mismatched_relationship
#
# Note: This script is optional. The relationship discrimination question is
# partially addressed by the FPR cutoff analysis (plots_proportion_exceeding_cutoffs.sh).
# Use this script if full continuous LR distributions are needed in addition
# to the threshold-based classification rates.

# Create logs directory if it doesn't exist
mkdir -p logs

# Load required modules
module load Rtidyverse

# Get input directory from command line argument
if [ -z "$1" ]; then
    echo "ERROR: Input directory required"
    echo "Usage: sbatch code/plots_mismatched_relationship.sh <input_dir> [output_subdir]"
    exit 1
fi

INPUT_DIR=$1

# Check if input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory does not exist: $INPUT_DIR"
    exit 1
fi

# Get optional output subdirectory
OUTPUT_SUBDIR=${2:-""}

echo "Input directory: $INPUT_DIR"
echo "Output directory: ${OUTPUT_SUBDIR:-"$INPUT_DIR/publication_figures (default)"}"

# Run the relationship discrimination plotting script
if [ -n "$OUTPUT_SUBDIR" ]; then
    Rscript code/plots_mismatched_relationship.R "$INPUT_DIR" "$OUTPUT_SUBDIR"
else
    Rscript code/plots_mismatched_relationship.R "$INPUT_DIR"
fi
