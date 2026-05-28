# ------------------------------------------------------------------------------
# Module 11: Ranking Test LR Calculator
# ------------------------------------------------------------------------------
#
# Purpose:
#   Takes assembled ranking databases from Module 10 and calculates combined
#   LRs for each focal-vs-database member pair across all loci sets, using
#   Module 4 (single locus LR) and Module 5 (combined LR).
#
#   Results are saved to output/focal_ranking_test/ and returned as a list
#   for direct input to Module 11 (ranking and outcome recorder).
#
# Dependencies:
#   - Module 4: Single Locus LR Calculator
#   - Module 5: Combined LR Calculator
#   - Module 10 output: assembled ranking databases
#
# ------------------------------------------------------------------------------

library(dplyr)
library(data.table)
library(purrr)

source("code/LR_kinship_utility_functions.R")
source("code/module1_allele_simulator.R")
source("code/module2_STR_profile_simulator.R")
source("code/module3_related_individual_simulator.R")
source("code/module4_single_locus_LR.R")
source("code/module5_combined_LR.R")


#' Calculate combined LRs for a single ranking search database
#'
#' Takes one assembled paired database (from Module 10) and calculates
#' combined LRs across all loci sets for the specified tested relationships
#' and population frequencies.
#'
#' @param paired_db Data frame, output from Module 10 assemble_ranking_database().
#'        Must have columns: batch_id, pair_id, individual_id, population, locus,
#'        focal_A1, focal_A2, ind2_A1, ind2_A2, known_relationship, is_true_relative
#' @param tested_relationships Character vector, relationships to test as H1.
#'        Default is c("parent_child", "full_siblings") per Tina's direction.
#' @param tested_populations Character vector, population frequencies to use.
#'        Default is "all" for the test run.
#' @param loci_sets Named list of loci vectors for combined LR calculation.
#' @param allele_frequency_data Data frame with allele frequency data.
#' @param kinship_coefficients Data frame with kinship coefficients.
#' @param output_dir Character, directory to save results.
#' @param save_results Logical, whether to save results to CSV (default TRUE).
#' @return Data frame of combined LRs with is_true_relative flag preserved,
#'         one row per pair per loci set per tested relationship x population.

calculate_ranking_lrs <- function(paired_db,
                                  tested_relationships = c("parent_child", "full_siblings"),
                                  tested_populations = "all",
                                  loci_sets,
                                  allele_frequency_data,
                                  kinship_coefficients,
                                  output_dir = "output/focal_ranking_test",
                                  save_results = TRUE) {

  # --- Input validation ---
  required_cols <- c("batch_id", "pair_id", "individual_id", "population", "locus",
                     "focal_A1", "focal_A2", "ind2_A1", "ind2_A2",
                     "known_relationship", "is_true_relative")
  missing_cols <- setdiff(required_cols, names(paired_db))
  if (length(missing_cols) > 0)
    stop(paste("paired_db missing columns:", paste(missing_cols, collapse = ", ")))

  if (!is.list(loci_sets) || is.null(names(loci_sets)))
    stop("loci_sets must be a named list")

  # Preserve is_true_relative and individual_id for later joining
  # Module 4 doesn't know about these so we keep them on the side
  pair_metadata <- paired_db |>
    select(pair_id, individual_id, known_relationship, is_true_relative) |>
    distinct()

  batch_id <- unique(paired_db$batch_id)
  cat(sprintf("Calculating LRs for batch: %s\n", batch_id))
  cat(sprintf("  Database size: %d pairs\n", length(unique(paired_db$pair_id))))
  cat(sprintf("  Tested relationships: %s\n", paste(tested_relationships, collapse = ", ")))
  cat(sprintf("  Tested populations: %s\n", paste(tested_populations, collapse = ", ")))
  cat(sprintf("  Loci sets: %s\n", paste(names(loci_sets), collapse = ", ")))

  # --- Calculate single locus LRs for each tested relationship ---
  # Module 4 takes one tested_relationship at a time
  single_locus_all <- map_dfr(tested_relationships, function(rel) {
    cat(sprintf("  Running single locus LR for tested relationship: %s\n", rel))
    calculate_single_locus_lr(
      pair_data            = paired_db,
      tested_relationship  = rel,
      tested_populations   = tested_populations,
      allele_frequency_data = allele_frequency_data,
      kinship_coefficients = kinship_coefficients
    )
  })

  # --- Calculate combined LRs across loci sets (Module 5) ---
  cat("  Calculating combined LRs across loci sets...\n")
  combined_lrs <- calculate_combined_lr(
    single_locus_results = single_locus_all,
    loci_sets            = loci_sets
  )

  # --- Rejoin metadata (individual_id, is_true_relative) ---
  combined_lrs <- combined_lrs |>
    left_join(pair_metadata, by = c("pair_id", "known_relationship"))

  # --- Save results ---
  if (save_results) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    filename  <- paste0("ranking_lrs_", batch_id, ".csv")
    file_path <- file.path(output_dir, filename)
    fwrite(combined_lrs, file_path)
    cat(sprintf("  Saved: %s\n", file_path))
  }

  return(combined_lrs)
}


#' Calculate combined LRs for all replicates
#'
#' Loops over a list of assembled ranking databases (from
#' assemble_ranking_replicates()) and calculates combined LRs for each,
#' returning a single combined data frame across all replicates.
#'
#' @param replicate_dbs Named list of paired data frames from Module 10
#'        assemble_ranking_replicates().
#' @param tested_relationships Character vector, relationships to test as H1.
#' @param tested_populations Character vector, population frequencies to use.
#' @param loci_sets Named list of loci vectors.
#' @param allele_frequency_data Data frame with allele frequency data.
#' @param kinship_coefficients Data frame with kinship coefficients.
#' @param output_dir Character, directory to save results.
#' @param save_combined Logical, whether to save the combined results CSV (default TRUE).
#' @return Data frame of combined LRs across all replicates.

calculate_ranking_lrs_replicates <- function(replicate_dbs,
                                             tested_relationships = c("parent_child", "full_siblings"),
                                             tested_populations = "all",
                                             loci_sets,
                                             allele_frequency_data,
                                             kinship_coefficients,
                                             output_dir = "output/focal_ranking_test",
                                             save_combined = TRUE) {

  n_replicates <- length(replicate_dbs)
  cat(sprintf("Calculating LRs for %d replicates...\n", n_replicates))

  all_lr_results <- imap_dfr(replicate_dbs, function(paired_db, focal_id) {
    cat(sprintf("\n--- Replicate focal_id: %s ---\n", focal_id))
    result <- calculate_ranking_lrs(
      paired_db             = paired_db,
      tested_relationships  = tested_relationships,
      tested_populations    = tested_populations,
      loci_sets             = loci_sets,
      allele_frequency_data = allele_frequency_data,
      kinship_coefficients  = kinship_coefficients,
      output_dir            = output_dir,
      save_results          = FALSE  # Save combined at the end instead
    )
    result$focal_id <- focal_id
    result
  })

  # --- Save combined results across all replicates ---
  if (save_combined) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    datetime_str <- format(Sys.time(), "%Y%m%d_%H%M%S")
    filename     <- paste0("ranking_lrs_all_replicates_", datetime_str, ".csv")
    file_path    <- file.path(output_dir, filename)
    fwrite(all_lr_results, file_path)
    cat(sprintf("\nAll replicates combined and saved: %s\n", file_path))
  }

  return(all_lr_results)
}


# ------------------------------------------------------------------------------
# Usage Examples (commented out)
# ------------------------------------------------------------------------------
#
# # --- Define loci sets (all 29 loci, subsets defined here for flexibility) ---
# loci_sets <- list(
#   core_13      = core_13_loci,
#   expanded_20  = expanded_20_loci,
#   autosomal_29 = all_29_loci
# )
#
# # --- Single replicate ---
# lr_results <- calculate_ranking_lrs(
#   paired_db             = paired_db,         # from Module 10
#   tested_relationships  = c("parent_child", "full_siblings"),
#   tested_populations    = "all",
#   loci_sets             = loci_sets,
#   allele_frequency_data = df_allelefreq,
#   kinship_coefficients  = kinship_matrix,
#   output_dir            = "output/focal_ranking_test"
# )
#
# # --- All replicates ---
# all_lr_results <- calculate_ranking_lrs_replicates(
#   replicate_dbs         = replicate_dbs,     # from Module 10
#   tested_relationships  = c("parent_child", "full_siblings"),
#   tested_populations    = "all",
#   loci_sets             = loci_sets,
#   allele_frequency_data = df_allelefreq,
#   kinship_coefficients  = kinship_matrix,
#   output_dir            = "output/focal_ranking_test"
# )
#
# # --- Then pass to Module 11 (ranking and outcome recorder) ---
# # ranking_results <- record_ranking_outcomes(all_lr_results, top_n = 200)
