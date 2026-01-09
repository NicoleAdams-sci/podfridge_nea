#!/bin/bash
#SBATCH --job-name=plot_cutoff
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=35G
#SBATCH --time=01:00:00
#SBATCH --output=logs/plot_cutoff_%j.out
##SBATCH --error=logs/plot_cutoff_%j.err

# Usage: sbatch code/plots_proportion_exceeding_cutoffs.sh <input_dir> [output_dir]

input_dir=$1
output_dir=$2

module load Rtidyverse
Rscript code/plots_proportion_exceeding_cutoffs.R $input_dir $output_dir
