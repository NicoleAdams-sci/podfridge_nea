#!/bin/bash
#SBATCH --job-name=sim_analysis
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=35G
#SBATCH --time=01:00:00
#SBATCH --output=logs/sim_analysis_%j.out
##SBATCH --error=logs/sim_analysis_%j.err

# Usage: sbatch code/simulation_analysis.sh <input_dir> [output_dir]

input_dir=$1
output_dir=${2:-$input_dir}  # Use input_dir as default

module load Rtidyverse
Rscript -e "rmarkdown::render('analysis/simulation_analysis.Rmd', \
  params = list(input_file = \"../output/$input_dir/combined_LR_all.rds\", \
                use_toy_data = FALSE, \
                output_dir = \"../output/$output_dir/analysis_results\"))"
