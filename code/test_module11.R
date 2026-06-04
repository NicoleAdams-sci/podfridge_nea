##### Module 11 Test — Single Replicate #####
# USAGE: Rscript code/test_module11.R <FOCAL_ID> <FOCAL_FAMILY_FILE> <FOCAL_POPULATION> <TESTED_POPULATIONS> <RELATIVE_TYPE>
#
# Examples:
#   Rscript code/test_module11.R 1 output/focal_ranking_test/Asian_focal_parent_child1_full_siblings1_20260604_123505.csv Asian all parent_child
#   Rscript code/test_module11.R 1 output/focal_ranking_test/Asian_focal_second_cousins1_20260604_170702.csv Asian all second_cousins

# ------------------------------------------------------------------------------
# 0. Parse arguments
# ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript test_module11.R <FOCAL_ID> <FOCAL_FAMILY_FILE> <FOCAL_POPULATION> <TESTED_POPULATIONS> <RELATIVE_TYPE>")
}

FOCAL_ID           <- args[1]
FOCAL_FAMILY_FILE  <- args[2]
FOCAL_POPULATION   <- args[3]
TESTED_POPULATIONS <- args[4]
RELATIVE_TYPE      <- args[5]
OUTPUT_DIR         <- "output/focal_ranking_test"

# Tested relationship = true relationship
TESTED_RELATIONSHIPS <- RELATIVE_TYPE

cat(sprintf("=== Module 11 Test — Replicate %s ===\n", FOCAL_ID))
cat("Started at:", format(Sys.time()), "\n\n")
cat(sprintf("  Focal ID             : %s\n", FOCAL_ID))
cat(sprintf("  Focal family file    : %s\n", basename(FOCAL_FAMILY_FILE)))
cat(sprintf("  Focal population     : %s\n", FOCAL_POPULATION))
cat(sprintf("  Tested populations   : %s\n", TESTED_POPULATIONS))
cat(sprintf("  Relative type        : %s\n\n", RELATIVE_TYPE))

# ------------------------------------------------------------------------------
# 1. Load packages and source modules
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(purrr)
})

source("code/LR_kinship_utility_functions.R")
source("code/module1_allele_simulator.R")
source("code/module2_STR_profile_simulator.R")
source("code/module3_related_individual_simulator.R")
source("code/module4_single_locus_LR.R")
source("code/module5_combined_LR.R")
source("code/module10_ranking_db_assembler.R")
source("code/module11_ranking_lr_calculator.R")

kinship_matrix <- fread("data/kinship_coefficients.csv")

# ------------------------------------------------------------------------------
# 2. Load focal family from explicit file path
# ------------------------------------------------------------------------------

cat("--- Loading focal family ---\n")

if (!file.exists(FOCAL_FAMILY_FILE))
  stop(paste("Focal family file not found:", FOCAL_FAMILY_FILE))

focal_family_data <- fread(FOCAL_FAMILY_FILE)

# Confirm relative type exists in file
available_rels <- unique(focal_family_data$relationship_to_focal)
if (!RELATIVE_TYPE %in% available_rels)
  stop(sprintf("Relative type '%s' not found in focal family file.\nAvailable: %s",
               RELATIVE_TYPE, paste(available_rels, collapse = ", ")))

# Confirm focal_id exists
available_ids <- unique(focal_family_data$focal_id[
  focal_family_data$relationship_to_focal == "self"
])
if (!FOCAL_ID %in% as.character(available_ids))
  stop(sprintf("focal_id %s not found. Available: %s",
               FOCAL_ID, paste(available_ids, collapse = ", ")))

cat(sprintf("Focal ID %s confirmed. Relative type '%s' confirmed.\n\n",
            FOCAL_ID, RELATIVE_TYPE))

# ------------------------------------------------------------------------------
# 3. Load unrelated pool
# ------------------------------------------------------------------------------

cat("--- Loading unrelated pool ---\n")

pool_files <- list.files("output/unrelated_pool", pattern = "\\.csv$", full.names = TRUE)
if (length(pool_files) == 0)
  stop("No unrelated pool files found. Run generate_unrelated_pool.sh first.")

unrelated_pool_data <- map_dfr(pool_files, fread)
n_unrelated <- nrow(unrelated_pool_data) / length(unique(unrelated_pool_data$locus))
cat(sprintf("Total unrelated individuals: %d\n\n", n_unrelated))

# ------------------------------------------------------------------------------
# 4. Assemble ranking database (module 10)
# ------------------------------------------------------------------------------

cat("--- Assembling ranking database (module 10) ---\n")

paired_db <- assemble_ranking_database(
  focal_family_data   = focal_family_data,
  unrelated_pool_data = unrelated_pool_data,
  focal_id            = FOCAL_ID,
  relative_type       = RELATIVE_TYPE,
  relative_index      = 1
)
cat(sprintf("Assembly complete: %d rows\n\n", nrow(paired_db)))

# ------------------------------------------------------------------------------
# 5. Calculate LRs (module 11)
# ------------------------------------------------------------------------------

cat("--- Calculating LRs (module 11) ---\n")

lr_results <- calculate_ranking_lrs(
  paired_db             = paired_db,
  tested_relationships  = TESTED_RELATIONSHIPS,
  tested_populations    = TESTED_POPULATIONS,
  loci_sets             = loci_lists,
  allele_frequency_data = df_allelefreq,
  kinship_coefficients  = kinship_matrix,
  output_dir            = OUTPUT_DIR,
  save_results          = FALSE
)

lr_results$focal_id      <- FOCAL_ID
lr_results$true_rel_type <- RELATIVE_TYPE

# ------------------------------------------------------------------------------
# 6. Save output for this replicate
# ------------------------------------------------------------------------------

cat("\n--- Saving results ---\n")

out_file <- file.path(OUTPUT_DIR,
                      sprintf("lr_focal_%s_%s.csv", FOCAL_ID, RELATIVE_TYPE))
fwrite(lr_results, out_file)
cat(sprintf("Saved: %s\n", out_file))
cat(sprintf("Rows: %d\n", nrow(lr_results)))

# ------------------------------------------------------------------------------
# 7. Quick sanity check
# ------------------------------------------------------------------------------

true_rel_lr <- lr_results |>
  filter(is_true_relative == TRUE, loci_set == "autosomal_29") |>
  pull(combined_LR)

cat(sprintf("\nSanity check — true relative LR (%s, autosomal_29): %.2e\n",
            RELATIVE_TYPE, as.numeric(true_rel_lr)))

if (length(true_rel_lr) == 1 && !is.na(as.numeric(true_rel_lr))) {
  cat("SUCCESS: Replicate", FOCAL_ID, "completed.\n")
} else {
  cat("WARNING: Unexpected true relative LR. Review output.\n")
}

cat("Completed at:", format(Sys.time()), "\n")
