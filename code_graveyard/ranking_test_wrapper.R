##### Focal Individual Ranking Test — R Wrapper #####
# USAGE: Rscript code/ranking_test_wrapper.R <N_FOCAL> <N_UNRELATED_PER_POP> <FOCAL_POP> <TESTED_POP> <TOP_N>
# Example: Rscript code/ranking_test_wrapper.R 10 250 Asian all 200

# ------------------------------------------------------------------------------
# 0. Parse arguments
# ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript ranking_test_wrapper.R <N_FOCAL> <N_UNRELATED_PER_POP> <FOCAL_POP> <TESTED_POP> <TOP_N>")
}

N_FOCAL_REPLICATES  <- as.integer(args[1])
N_UNRELATED_PER_POP <- as.integer(args[2])
FOCAL_POPULATION    <- args[3]
TESTED_POPULATIONS  <- args[4]
TOP_N               <- as.integer(args[5])

UNRELATED_POOL_DIR  <- "output/unrelated_pool"
RANKING_TEST_DIR    <- "output/focal_ranking_test"
RELATIVE_TYPES      <- c("parent_child", "full_siblings")

cat("Started at:", format(Sys.time()), "\n")

# ------------------------------------------------------------------------------
# 1. Load packages and source modules
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(purrr)
})

# Utility functions — also defines df_allelefreq, loci_lists, loci_list
source("code/LR_kinship_utility_functions.R")

# Simulation modules
source("code/module1_allele_simulator.R")
source("code/module2_STR_profile_simulator.R")
source("code/module3_related_individual_simulator.R")
source("code/module4_single_locus_LR.R")
source("code/module5_combined_LR.R")
source("code/module7_single_pop_focal_generator.R")
source("code/module8_unrelated_pool_generator.R")

# Ranking pipeline modules
source("code/module10_ranking_db_assembler.R")
source("code/module11_ranking_lr_calculator.R")
source("code/module12_ranking_outcome_recorder.R")

# Load kinship coefficients
# Note: must be named kinship_matrix — module 3 references this as a global variable
kinship_matrix <- fread("data/kinship_coefficients.csv")

# ------------------------------------------------------------------------------
# STEP 1: Generate unrelated pool (once, reused across all replicates)
# ------------------------------------------------------------------------------

cat("\n--- STEP 1: Generate unrelated pool ---\n")

unrelated_results <- generate_multiple_pop_unrelated(
  populations           = c("AfAm", "Cauc", "Hispanic", "Asian"),
  n_unrelated_per_pop   = N_UNRELATED_PER_POP,
  loci_list             = loci_list,
  allele_frequency_data = df_allelefreq,
  output_dir            = UNRELATED_POOL_DIR,
  use_single_datetime   = TRUE
)

# Load and combine all population files into one pool
unrelated_pool_data <- map_dfr(unrelated_results$file_path, fread)
cat(sprintf("Total unrelated individuals loaded: %d\n",
            length(unique(unrelated_pool_data$individual_id))))

# ------------------------------------------------------------------------------
# STEP 2: Generate focal individuals and relatives (Module 7)
# ------------------------------------------------------------------------------

cat("\n--- STEP 2: Generate focal individuals and relatives ---\n")

focal_result <- generate_single_pop_focal(
  population            = FOCAL_POPULATION,
  n_focal               = N_FOCAL_REPLICATES,
  relationship_counts   = list(parent_child = 1, full_siblings = 1),
  loci_list             = loci_list,
  allele_frequency_data = df_allelefreq,
  kinship_coefficients  = kinship_matrix,
  output_dir            = RANKING_TEST_DIR
)

focal_family_data <- focal_result$data
cat(sprintf("Focal individuals generated: %d\n",
            length(unique(focal_family_data$focal_id[
              focal_family_data$relationship_to_focal == "self"
            ]))))

# ------------------------------------------------------------------------------
# STEP 3: Assemble ranking databases (Module 10)
# ------------------------------------------------------------------------------

cat("\n--- STEP 3: Assemble ranking databases ---\n")

replicate_dbs <- map(RELATIVE_TYPES, function(rel_type) {
  cat(sprintf("\nAssembling databases for: %s\n", rel_type))
  assemble_ranking_replicates(
    focal_family_data   = focal_family_data,
    unrelated_pool_data = unrelated_pool_data,
    relative_type       = rel_type,
    relative_index      = 1
  )
}) |> set_names(RELATIVE_TYPES)

# ------------------------------------------------------------------------------
# STEP 4: Calculate combined LRs (Module 11)
# ------------------------------------------------------------------------------

cat("\n--- STEP 4: Calculate combined LRs ---\n")

all_lr_results <- map_dfr(RELATIVE_TYPES, function(rel_type) {
  cat(sprintf("\nCalculating LRs for: %s\n", rel_type))
  result <- calculate_ranking_lrs_replicates(
    replicate_dbs         = replicate_dbs[[rel_type]],
    tested_relationships  = c("parent_child", "full_siblings"),
    tested_populations    = TESTED_POPULATIONS,
    loci_sets             = loci_lists,
    allele_frequency_data = df_allelefreq,
    kinship_coefficients  = kinship_matrix,
    output_dir            = RANKING_TEST_DIR,
    save_combined         = TRUE
  )
  result$true_relative_type <- rel_type
  result
})

cat(sprintf("\nTotal LR rows calculated: %s\n",
            format(nrow(all_lr_results), big.mark = ",")))

# ------------------------------------------------------------------------------
# STEP 5: Rank and record outcomes (Module 12)
# ------------------------------------------------------------------------------

cat("\n--- STEP 5: Rank and record outcomes ---\n")

ranking_results <- record_ranking_outcomes(
  lr_results   = all_lr_results,
  top_n        = TOP_N,
  output_dir   = RANKING_TEST_DIR,
  save_results = TRUE
)

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------

cat("\n=============================================================\n")
cat("  RANKING TEST COMPLETE\n")
cat("=============================================================\n")
cat(sprintf("  Outputs saved to: %s\n\n", RANKING_TEST_DIR))
cat("  Results summary:\n")
print(ranking_results$summary |>
        select(loci_set, tested_relationship, known_relationship,
               n_replicates, prop_in_top_n, median_rank))

cat("\nCompleted at:", format(Sys.time()), "\n")
cat("SUCCESS: Ranking test completed\n")
