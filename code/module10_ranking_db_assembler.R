# ------------------------------------------------------------------------------
# Module 10: Ranking Database Assembler
# ------------------------------------------------------------------------------
#
# Purpose:
#   Takes a focal individual + one relative (from Module 7 output) and an
#   unrelated pool (from Module 8 output), and assembles them into the paired
#   format required by Module 4 (calculate_single_locus_lr).
#
#   The focal individual's genotype is held constant and paired against every
#   individual in the database (unrelated pool + one true relative), producing
#   one row per focal-database member pair per locus.
#
# Dependencies:
#   - Module 7 output: focal individual + relatives
#   - Module 8 output: unrelated pool
#
# Output format matches Module 4 input requirements:
#   batch_id, pair_id, population, locus,
#   focal_A1, focal_A2, ind2_A1, ind2_A2, known_relationship
#
# ------------------------------------------------------------------------------

library(dplyr)
library(data.table)
library(purrr)

#' Assemble a ranking search database for one focal individual and one relative type
#'
#' Takes one focal individual's genotype, one of their simulated relatives, and
#' the full unrelated pool, then pairs the focal against every database member
#' (relative + unrelated) in the format Module 4 expects.
#'
#' @param focal_family_data Data frame, output from Module 7 for one focal individual.
#'        Must contain rows where relationship_to_focal == "self" (the focal)
#'        and at least one relative row.
#' @param unrelated_pool_data Data frame, output from Module 8.
#'        Contains all unrelated individuals across populations.
#' @param focal_id Character, the focal_id to use from focal_family_data
#'        (e.g. "001"). If NULL, uses the first focal individual found.
#' @param relative_type Character, which relative to insert into the database.
#'        Must match a value in relationship_to_focal column of focal_family_data.
#'        (e.g. "parent_child", "full_siblings")
#' @param relative_index Integer, which replicate of that relative type to use
#'        if multiple were simulated (default = 1).
#' @param batch_id Character, batch identifier for this search replicate.
#'        Defaults to current datetime.
#' @return Data frame in Module 4 format, one row per locus per database member

assemble_ranking_database <- function(focal_family_data,
                                      unrelated_pool_data,
                                      focal_id = NULL,
                                      relative_type,
                                      relative_index = 1,
                                      batch_id = NULL) {

  # --- Input validation ---
  required_focal_cols <- c("family_id", "focal_id", "individual_id",
                            "relationship_to_focal", "population", "locus", "A1", "A2")
  required_unrel_cols <- c("individual_id", "relationship_to_focal",
                            "population", "locus", "A1", "A2")

  missing_focal <- setdiff(required_focal_cols, names(focal_family_data))
  missing_unrel  <- setdiff(required_unrel_cols, names(unrelated_pool_data))

  if (length(missing_focal) > 0)
    stop(paste("focal_family_data missing columns:", paste(missing_focal, collapse = ", ")))
  if (length(missing_unrel) > 0)
    stop(paste("unrelated_pool_data missing columns:", paste(missing_unrel, collapse = ", ")))

  # --- Generate batch_id if not provided ---
  if (is.null(batch_id)) {
    batch_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
  }

  # --- Select focal individual ---
  if (is.null(focal_id)) {
    focal_id <- focal_family_data$focal_id[focal_family_data$relationship_to_focal == "self"][1]
    message(paste("No focal_id provided. Using first found:", focal_id))
  }

  focal_genotype <- focal_family_data |>
    filter(focal_id == !!focal_id, relationship_to_focal == "self") |>
    select(locus, focal_A1 = A1, focal_A2 = A2, population)

  if (nrow(focal_genotype) == 0)
    stop(paste("No focal individual found with focal_id:", focal_id))

  focal_population <- unique(focal_genotype$population)

  # --- Select the one relative to insert ---
  # relative individual_id pattern: e.g. "parentchild1_001", "fullsiblings2_001"
  available_relatives <- focal_family_data |>
    filter(focal_id == !!focal_id,
           relationship_to_focal == relative_type)

  if (nrow(available_relatives) == 0)
    stop(paste("No relative of type", relative_type, "found for focal_id:", focal_id))

  # Get the relative_index-th unique individual of this type
  relative_ids <- unique(available_relatives$individual_id)
  if (relative_index > length(relative_ids))
    stop(paste("relative_index", relative_index, "exceeds number of available relatives:",
               length(relative_ids)))

  selected_relative_id <- relative_ids[relative_index]

  relative_genotype <- available_relatives |>
    filter(individual_id == selected_relative_id) |>
    select(locus, ind2_A1 = A1, ind2_A2 = A2) |>
    mutate(ind2_A1 = as.character(ind2_A1),
           ind2_A2 = as.character(ind2_A2),
           individual_id      = selected_relative_id,
           known_relationship = relative_type,
           is_true_relative   = TRUE)

  # --- Format unrelated pool ---
  # Prefix individual_id with population to avoid ID collisions when
  # module 8 generates the same IDs (unrel_001 etc.) across populations
  unrelated_genotype <- unrelated_pool_data |>
    mutate(individual_id = paste0(population, "_", individual_id)) |>
    select(locus, ind2_A1 = A1, ind2_A2 = A2, individual_id) |>
    mutate(ind2_A1 = as.character(ind2_A1),
           ind2_A2 = as.character(ind2_A2),
           known_relationship = "unrelated",
           is_true_relative   = FALSE)

  # --- Combine relative + unrelated into one database ---
  database <- bind_rows(relative_genotype, unrelated_genotype)

  n_database  <- length(unique(database$individual_id))
  n_unrelated <- n_database - 1  # subtract the one true relative
  cat(sprintf("Database assembled: %d individuals (%d unrelated + 1 %s)\n",
              n_database, n_unrelated, relative_type))

  # --- Join focal genotype against every database member ---
  # One row per locus per database member
  paired_data <- database |>
    left_join(focal_genotype, by = "locus") |>
    mutate(
      batch_id   = batch_id,
      pair_id    = sprintf("%04d", as.integer(factor(individual_id))),
      # population is the focal individual's true population
      population = focal_population
    ) |>
    select(batch_id, pair_id, individual_id, population, locus,
           focal_A1, focal_A2, ind2_A1, ind2_A2,
           known_relationship, is_true_relative)

  # Verify all loci present for all pairs
  n_loci      <- length(unique(paired_data$locus))
  n_pairs_out <- length(unique(paired_data$pair_id))
  cat(sprintf("Paired data: %d pairs x %d loci = %d rows\n",
              n_pairs_out, n_loci, nrow(paired_data)))

  return(paired_data)
}


#' Assemble ranking databases for multiple replicates
#'
#' Loops over multiple focal individuals (replicates) and assembles a search
#' database for each, returning a list of paired data frames ready for Module 4.
#'
#' @param focal_family_data Data frame, full Module 7 output (multiple focal individuals)
#' @param unrelated_pool_data Data frame, Module 8 output (fixed unrelated pool)
#' @param relative_type Character, relationship type to test (e.g. "parent_child")
#' @param relative_index Integer, which replicate of the relative to use (default = 1)
#' @return Named list of paired data frames, one per focal individual (replicate)

assemble_ranking_replicates <- function(focal_family_data,
                                        unrelated_pool_data,
                                        relative_type,
                                        relative_index = 1) {

  # Get all focal IDs
  focal_ids <- unique(focal_family_data$focal_id[
    focal_family_data$relationship_to_focal == "self"
  ])

  cat(sprintf("Assembling ranking databases for %d replicates (%s)...\n",
              length(focal_ids), relative_type))

  # Shared batch_id across all replicates so they're linkable
  shared_batch_id <- format(Sys.time(), "%Y%m%d_%H%M%S")

  replicate_databases <- map(focal_ids, function(fid) {
    assemble_ranking_database(
      focal_family_data  = focal_family_data,
      unrelated_pool_data = unrelated_pool_data,
      focal_id           = fid,
      relative_type      = relative_type,
      relative_index     = relative_index,
      batch_id           = paste0(shared_batch_id, "_focal", fid)
    )
  })

  names(replicate_databases) <- focal_ids

  cat(sprintf("Done. %d replicate databases assembled.\n", length(replicate_databases)))

  return(replicate_databases)
}


# ------------------------------------------------------------------------------
# Usage Examples (commented out)
# ------------------------------------------------------------------------------
#
# # --- Load data ---
# # Unrelated pool lives in its own directory (generated once, reused across all tests)
# unrelated_pool_data <- fread("output/unrelated_pool/unrelated_pool_20250528.csv")
#
# # Focal families live in the ranking test directory
# focal_family_data <- fread("output/focal_ranking_test/Asian_focal_parentchild1_fullsiblings1_20250528.csv")
#
# # --- Single replicate (one focal individual, testing parent_child) ---
# paired_db <- assemble_ranking_database(
#   focal_family_data   = focal_family_data,
#   unrelated_pool_data = unrelated_pool_data,
#   focal_id            = "001",
#   relative_type       = "parent_child"
# )
#
# # --- All replicates for parent_child ---
# replicate_dbs <- assemble_ranking_replicates(
#   focal_family_data   = focal_family_data,
#   unrelated_pool_data = unrelated_pool_data,
#   relative_type       = "parent_child"
# )
#
# # --- Then pass to Module 11 (LR calculator) ---
# # lr_results <- calculate_ranking_lrs(replicate_dbs, ...)
