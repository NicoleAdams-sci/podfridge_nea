# ------------------------------------------------------------------------------
# Module 6: Pair Batch Generator
# ------------------------------------------------------------------------------
#
# Dependencies:
#   - Module 1: Allele Simulator (via Module 3)
#   - Module 2: STR Profile Simulator (via Module 3)
#   - Module 3: Related Individual Simulator (direct)
#   - LR_kinship_utility_functions.R (via Module 3)
#
# Direct function calls:
#   - simulate_multiple_pairs() from Module 3
#
# ------------------------------------------------------------------------------

library(dplyr)
library(data.table)

# Source required modules
source("code/LR_kinship_utility_functions.R")
source("code/module1_allele_simulator.R")
source("code/module2_STR_profile_simulator.R")
source("code/module3_related_individual_simulator.R")

#' Generate pairs for a single population-relationship combination and save to CSV
#'
#' Creates pairs for one specific population and relationship combination,
#' then saves the results to a CSV file with standardized naming
#'
#' @param population Character, population for simulation (e.g., "AfAm", "Cauc")
#' @param relationship Character, relationship type (e.g., "full_siblings", "parent_child")
#' @param n_pairs Integer, number of pairs to simulate
#' @param loci_list Character vector, loci to simulate (defaults to all available)
#' @param allele_frequency_data Data frame with allele frequency data
#' @param kinship_coefficients Data frame with kinship coefficients
#' @param output_dir Character, directory to save files (defaults to "output")
#' @param custom_datetime Character, custom datetime string (defaults to current time)
#' @return List with 'data' (the pairs data frame) and 'file_path' (path to saved CSV)
generate_pair_batch <- function(population, 
                                         relationship, 
                                         n_pairs,
                                         loci_list = NULL,
                                         allele_frequency_data,
                                         kinship_coefficients,
                                         output_dir = "output",
                                         custom_datetime = NULL) {
  
  # Validate inputs
  if (!population %in% unique(allele_frequency_data$population)) {
    stop(paste("Population", population, "not found in allele_frequency_data"))
  }
  
  if (!relationship %in% kinship_coefficients$relationship_type) {
    stop(paste("Relationship", relationship, "not found in kinship_coefficients"))
  }
  
  if (n_pairs < 1) {
    stop("n_pairs must be at least 1")
  }
  
  # Use all available loci if not specified
  if (is.null(loci_list)) {
    loci_list <- unique(allele_frequency_data$marker)
  }
  
  # Generate datetime string
  if (is.null(custom_datetime)) {
    datetime_str <- format(Sys.time(), "%Y%m%d_%H%M%S")
  } else {
    datetime_str <- custom_datetime
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Generate filename
  filename <- paste0(population, "_", relationship, "_", datetime_str, ".csv")
  file_path <- file.path(output_dir, filename)
  
  # Generate pairs using Module 3
  cat("Generating", n_pairs, "pairs for", population, "-", relationship, "...\n")
  
  pairs_data <- simulate_multiple_pairs(
    n_pairs = n_pairs,
    loci_list = loci_list,
    population = population,
    known_relationship = relationship,
    allele_frequency_data = allele_frequency_data,
    batch_id = datetime_str
  )
  
  # Save to CSV
  cat("Saving results to:", file_path, "\n")
  fwrite(pairs_data, file_path)
  
  # Print summary
  cat("Successfully generated:\n")
  cat("  - Pairs:", n_pairs, "\n")
  cat("  - Loci:", length(loci_list), "\n")
  cat("  - Total rows:", nrow(pairs_data), "\n")
  cat("  - File:", filename, "\n")
  
  return(list(
    data = pairs_data,
    file_path = file_path
  ))
}

#' Generate pairs for multiple population-relationship combinations
#'
#' Convenience function to generate multiple single combo files at once
#'
#' @param combinations Data frame with columns: population, relationship, n_pairs
#' @param loci_list Character vector, loci to simulate
#' @param allele_frequency_data Data frame with allele frequency data
#' @param kinship_coefficients Data frame with kinship coefficients
#' @param output_dir Character, directory to save files
#' @param use_single_datetime Logical, whether to use same datetime for all files
#' @return Data frame with combination details and file paths
generate_multiple_pair_batches <- function(combinations,
                                           loci_list = NULL,
                                           allele_frequency_data,
                                           kinship_coefficients,
                                           output_dir = "output",
                                           use_single_datetime = TRUE) {
  
  # Validate combinations input
  required_cols <- c("population", "relationship", "n_pairs")
  missing_cols <- setdiff(required_cols, names(combinations))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns in combinations:", paste(missing_cols, collapse = ", ")))
  }
  
  # Generate single datetime if requested
  if (use_single_datetime) {
    shared_datetime <- format(Sys.time(), "%Y%m%d_%H%M%S")
  } else {
    shared_datetime <- NULL
  }
  
  # Generate pairs for each combination
  results <- combinations %>%
    rowwise() %>%
    mutate(
      result = list(generate_pair_batch(
        population = population,
        relationship = relationship,
        n_pairs = n_pairs,
        loci_list = loci_list,
        allele_frequency_data = allele_frequency_data,
        kinship_coefficients = kinship_coefficients,
        output_dir = output_dir,
        custom_datetime = shared_datetime
      )),
      file_path = map_chr(result, ~ .x$file_path),
      filename = basename(file_path),
      datetime_used = ifelse(is.null(shared_datetime), 
                            format(Sys.time(), "%Y%m%d_%H%M%S"), 
                            shared_datetime)
    ) %>%
    select(-result) %>%
    ungroup()
  
  cat("\nGenerated", nrow(results), "combination files.\n")
  
  return(as.data.frame(results))
}

#' Create standard population-relationship combinations
#'
#' Helper function to create a combinations data frame for common scenarios
#'
#' @param populations Character vector of populations
#' @param relationships Character vector of relationships  
#' @param n_pairs_per_combo Integer, number of pairs for each combination
#' @return Data frame with all combinations
create_pair_combinations <- function(populations, relationships, n_pairs_per_combo) {
  
  combinations <- expand.grid(
    population = populations,
    relationship = relationships,
    stringsAsFactors = FALSE
  ) %>%
    mutate(n_pairs = n_pairs_per_combo)
  
  return(combinations)
}

# ------------------------------------------------------------------------------
# Usage Examples (commented out)
# ------------------------------------------------------------------------------
# # Load required data
# df_allelefreq <- fread("data/df_allelefreq_combined.csv")
# df_allelefreq <- df_allelefreq[population != "all"]
# df_allelefreq$frequency <- ifelse(df_allelefreq$frequency == 0, FALLBACK_FREQ, df_allelefreq$frequency)
# kinship_matrix <- fread("data/kinship_coefficients.csv")
# 
# # Define loci
# core_loci <- fread("data/core_CODIS_loci.csv")
# loci_list <- core_loci[core_loci$core_13 == 1]$locus
# 
# # Single combination
# result <- generate_pair_batch(
#   population = "AfAm",
#   relationship = "full_siblings", 
#   n_pairs = 100,
#   loci_list = loci_list,
#   allele_frequency_data = df_allelefreq,
#   kinship_coefficients = kinship_matrix
# )
# 
# # Access the data and file path
# pairs_data <- result$data
# file_path <- result$file_path
# 
# # Multiple combinations
# combinations <- create_pair_combinations(
#   populations = c("AfAm", "Cauc", "Hispanic", "Asian"),
#   relationships = c("parent_child", "full_siblings", "unrelated"),
#   n_pairs_per_combo = 50
# )
# 
# results <- generate_multiple_pair_batches(
#   combinations = combinations,
#   loci_list = loci_list,
#   allele_frequency_data = df_allelefreq,
#   kinship_coefficients = kinship_matrix
# )