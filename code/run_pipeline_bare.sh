#!/bin/bash
################################################################################
# PODFRIDGE Pipeline - Bare Commands Only
# 
# Copy and paste these commands one section at a time.
# Wait for each section to complete before running the next.
################################################################################

# ============================================================================
# STEP 1: SIMULATE PAIRS
# ============================================================================
sbatch code/sim_pairs.sh
sbatch code/sim_pairs_unrelated.sh
# Wait for: squeue -u $USER shows no jobs
# Verify: ls output/pairs_*.csv | wc -l  (should be 1000)

# ============================================================================
# STEP 2: CALCULATE SINGLE-LOCUS LRs
# ============================================================================
bash code/lr_submission.sh
# Wait for: squeue -u $USER shows no jobs
# Verify: ls output/LR/LR_*.csv | wc -l  (should be 1000)

# ============================================================================
# STEP 3: CALCULATE COMBINED LRs
# ============================================================================
bash code/combined_lr_submission.sh
# Wait for: squeue -u $USER shows no jobs
# Verify: ls output/combined_LR/combined_LR_*.csv | wc -l  (should be 1000)

# ============================================================================
# STEP 4: ANALYZE RESULTS
# ============================================================================
sbatch code/analyze_lr_outputs.sh output/lr_analysis_$(date +%Y%m%d)
# Wait for: squeue -u $USER shows no jobs
# Verify: ls output/lr_analysis_*/  (check directory created)

# ============================================================================
# STEP 5: GENERATE PLOTS
# ============================================================================
# Replace YYYYMMDD with your analysis date from Step 4

module load Rtidyverse

# Matched plots (run directly - fast)
Rscript code/plots_matched.R output/lr_analysis_YYYYMMDD lr_analysis_YYYYMMDD/plots_matched

# Mismatched plots (submit as job)
sbatch code/plots_mismatched.sh output/lr_analysis_YYYYMMDD lr_analysis_YYYYMMDD/plots_mismatched

# Threshold plots (submit as job)
sbatch code/plots_proportion_exceeding_cutoffs.sh lr_analysis_YYYYMMDD lr_analysis_YYYYMMDD/plots_exceeding_cutoffs

# Wait for: squeue -u $USER shows no jobs

# ============================================================================
# STEP 6: GENERATE STATISTICAL REPORT
# ============================================================================
sbatch code/simulation_analysis.sh lr_analysis_YYYYMMDD
# Wait for: squeue -u $USER shows no jobs
# Report at: output/lr_analysis_YYYYMMDD/analysis_results/simulation_analysis_YYYYMMDD.html

# ============================================================================
# DONE!
# ============================================================================
