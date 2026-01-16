#!/bin/bash
#SBATCH --job-name=sim_analysis
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=21G
#SBATCH --time=00:30:00
#SBATCH --output=logs/sim_analysis_%j.out
##SBATCH --error=logs/sim_analysis_%j.err

# Usage: sbatch code/simulation_analysis.sh <input_dir> [output_dir]

input_dir=$1
output_dir=${2:-$input_dir}  # Use input_dir as default

# Generate date string
date_str=$(date +%Y%m%d)

module load Rtidyverse
Rscript -e "rmarkdown::render('analysis/simulation_analysis.Rmd', \
  output_file = paste0('simulation_analysis_', '$date_str', '.html'), \
  params = list(input_file = \"../output/$input_dir/combined_LR_all.rds\", \
                use_toy_data = FALSE, \
                output_dir = \"../output/$output_dir/analysis_results\"))"
