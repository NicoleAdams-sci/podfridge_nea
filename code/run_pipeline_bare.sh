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
# STEP 2.5: ANALYZE PER-LOCUS INFLATION  (optional; runs in parallel w/ 3-6)
# ============================================================================
# Reads output/LR/ directly - no dependency on combined LRs or later steps.
sbatch code/analyze_locus_inflation.sh output/lr_analysis_YYYYMMDD/locus_inflation
# Verify: ls output/lr_analysis_YYYYMMDD/locus_inflation/*.csv

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
# STEP 4.5: PREPARE INTERMEDIATE CSVs FOR PLOTTING
# ============================================================================
sbatch code/prepare_combined_lr_intermediates.sh output/lr_analysis_YYYYMMDD
# Wait for: squeue -u $USER shows no jobs
# Verify: ls output/lr_analysis_YYYYMMDD/proportions_with_classification.csv

# ============================================================================
# STEP 5: GENERATE PUBLICATION PLOTS
# ============================================================================
# Replace YYYYMMDD with your analysis date from Step 4.
# All plotting scripts run interactively (Rscript), not via sbatch.

module load Rtidyverse

# Matched scenarios
Rscript code/plots_matched_publication.R output/lr_analysis_YYYYMMDD \
    output/lr_analysis_YYYYMMDD/plots_matched

# Population mismatch (needs Step 4.5)
Rscript code/plots_mismatched_population.R output/lr_analysis_YYYYMMDD \
    output/lr_analysis_YYYYMMDD/plots_mismatched_population

# Relationship discrimination (loads the full .rds - ~42 GB)
Rscript code/plots_mismatched_relationship.R output/lr_analysis_YYYYMMDD \
    output/lr_analysis_YYYYMMDD/plots_mismatched_relationship

# Classification / FPR thresholds (needs Step 4.5)
Rscript code/plots_cutoffs_publication.R output/lr_analysis_YYYYMMDD \
    output/lr_analysis_YYYYMMDD/plots_cutoffs

# Per-locus inflation (needs Step 2.5)
Rscript code/plots_locus_inflation.R output/lr_analysis_YYYYMMDD/locus_inflation \
    output/lr_analysis_YYYYMMDD/locus_inflation/plots

# ============================================================================
# STEP 6: RUN STATISTICAL TESTS
# ============================================================================
sbatch code/run_statistical_tests.sh output/lr_analysis_YYYYMMDD
# Wait for: squeue -u $USER shows no jobs
# Verify: ls output/lr_analysis_YYYYMMDD/stats/*.csv

# ============================================================================
# STEP 7: GENERATE STATISTICAL REPORT (optional)
# ============================================================================
sbatch code/simulation_analysis.sh lr_analysis_YYYYMMDD
# Wait for: squeue -u $USER shows no jobs
# Report at: output/lr_analysis_YYYYMMDD/analysis_results/simulation_analysis_YYYYMMDD.html

# ============================================================================
# FOCAL RANKING TEST (separate workflow)
# ============================================================================
# See README_focal_ranking_test.md for the full per-database-size workflow.

# ============================================================================
# DONE!
# ============================================================================
