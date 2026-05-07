#!/bin/bash
#SBATCH --job-name=locus_inflation
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/locus_inflation_%j.out
#SBATCH --error=logs/locus_inflation_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=nicolead@umich.edu

# SLURM wrapper for analyze_locus_inflation.R
#
# Identifies which loci drive LR inflation under mismatched population
# allele frequencies, using single-locus LR files from output/LR/.
#
# Usage:
#   sbatch code/analyze_locus_inflation.sh [output_dir]
#   sbatch code/analyze_locus_inflation.sh output/locus_inflation_20260130
#
# SLURM_CPUS_PER_TASK is passed automatically to the R script, which uses
# it to set the number of parallel workers via mclapply.

mkdir -p logs

module load Rtidyverse

echo "================================================"
echo "Locus Inflation Analysis"
echo "Job ID:    ${SLURM_JOB_ID}"
echo "Cores:     ${SLURM_CPUS_PER_TASK}"
echo "Started:   $(date)"
echo "================================================"

# Get optional output directory from command line
OUTPUT_DIR=${1:-""}

# Count input files for verification
LR_COUNT=$(ls output/LR/LR_*.csv 2>/dev/null | wc -l)
echo "Found ${LR_COUNT} single-locus LR files in output/LR/"
echo ""

# Run the analysis script
echo "Running analyze_locus_inflation.R..."
echo "================================================"

if [ -z "${OUTPUT_DIR}" ]; then
  Rscript code/analyze_locus_inflation.R
  EXITCODE=$?
else
  Rscript code/analyze_locus_inflation.R "${OUTPUT_DIR}"
  EXITCODE=$?
fi

echo "================================================"

if [ ${EXITCODE} -eq 0 ]; then
  echo "SUCCESS: Locus inflation analysis completed"

  # Find the output directory if not specified
  if [ -z "${OUTPUT_DIR}" ]; then
    OUTPUT_DIR=$(ls -dt output/locus_inflation_* 2>/dev/null | head -1)
  fi

  if [ -d "${OUTPUT_DIR}" ]; then
    echo "Output saved to: ${OUTPUT_DIR}"
    echo ""
    echo "Output files:"
    for f in locus_inflation_summary.csv locus_heterozygosity_summary.csv; do
      if [ -f "${OUTPUT_DIR}/$f" ]; then
        size=$(du -h "${OUTPUT_DIR}/$f" | cut -f1)
        echo "  - $f ($size)"
      fi
    done
  fi
else
  echo "ERROR: Analysis failed with exit code ${EXITCODE}"
  echo "Check logs:"
  echo "  logs/locus_inflation_${SLURM_JOB_ID}.out"
  echo "  logs/locus_inflation_${SLURM_JOB_ID}.err"
  exit ${EXITCODE}
fi

echo ""
echo "================================================"
echo "Completed: $(date)"
echo "Runtime:   ${SECONDS} seconds"
echo "================================================"
