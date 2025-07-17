# ------------------------------------------------------------------------------
# Module 5: Single Locus Likelihood Ratio Calculator
# ------------------------------------------------------------------------------

library(dplyr)
library(tibble)
library(purrr)

# Source the utility functions for kinship calculations
source("code/LR_kinship_utility_functions.R")

#' Calculate single-locus likelihood ratios for individual pairs
#'
#' This function takes genotype data for pairs of individuals and calculates
#' likelihood ratios for specified relationships and populations at each locus.
#'
#' @param pair_data Data frame with columns: pair_id, population, locus, 
#'                  focal_A1, focal_A2, ind2_A1, ind2_A2, known_relationship
#' @param tested_relationship Character, relationship to test (e.g., "full_siblings")
#' @param tested_populations Character vector, populations to test LRs against
#' @param allele_frequency_data Data frame with population, marker, allele, frequency
#' @param kinship_coefficients Data frame with relationship_type, k0, k1, k2
#' @return Data frame with LR calculations for each locus-population combination
calculate_single_locus_lr <- function(pair_data, 
                                      tested_relationship, 
                                      tested_populations,
                                      allele_frequency_data, 
                                      kinship_coefficients) {
  
  # Validate inputs
  required_cols <- c("batch_id", "pair_id", "population", "locus", "focal_A1", "focal_A2", 
                     "ind2_A1", "ind2_A2", "known_relationship")
  missing_cols <- setdiff(required_cols, names(pair_data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  if (!tested_relationship %in% kinship_coefficients$relationship_type) {
    stop(paste("tested_relationship", tested_relationship, 
               "not found in kinship_coefficients"))
  }
  
  # Standardize genotype column names for kinship_calculation function
  standardized_data <- pair_data |>
    mutate(
      ind1_allele1 = focal_A1,
      ind1_allele2 = focal_A2,
      ind2_allele1 = ind2_A1,
      ind2_allele2 = ind2_A2
    )
  
  # Calculate LRs for each row (locus) using the utility function
  results_list <- vector("list", nrow(standardized_data))
  
  for (i in 1:nrow(standardized_data)) {
    row <- standardized_data[i, ]
    
    # Ensure alleles are character vectors (not factors or lists)
    row$ind1_allele1 <- as.character(row$ind1_allele1)
    row$ind1_allele2 <- as.character(row$ind1_allele2)
    row$ind2_allele1 <- as.character(row$ind2_allele1)
    row$ind2_allele2 <- as.character(row$ind2_allele2)
    
    # Convert row to list format expected by kinship_calculation
    row_list <- as.list(row)
    
    # Use the kinship_calculation function from utility functions
    lr_results <- kinship_calculation(
      row = row_list,
      allele_frequency_data = allele_frequency_data,
      kinship_matrix = kinship_coefficients,
      tested_relationship = tested_relationship,
      tested_populations = tested_populations
    )
    
    results_list[[i]] <- lr_results
  }
  
  # Combine all results
  results <- bind_rows(results_list)
  
  # Clean up and reorganize output
  final_results <- results |>
    select(batch_id, pair_id, population, known_relationship, locus, 
           ind1_allele1, ind1_allele2, ind2_allele1, ind2_allele2,
           shared_alleles, genotype_match, tested_relationship, 
           tested_population, LR) |>
    rename(
      focal_A1 = ind1_allele1,
      focal_A2 = ind1_allele2,
      ind2_A1 = ind2_allele1,
      ind2_A2 = ind2_allele2
    )
  
  return(final_results)
}

# Note: Multiple relationships testing will be handled in a separate module if needed

# Note: Combined LR calculation across loci will be handled in a separate module

# Note: This function currently only processes one tested relationship at a time
# For batch processing multiple relationships, use separate calls or a future module

# ------------------------------------------------------------------------------
# Usage Example (commented out)
# ------------------------------------------------------------------------------
# # Load required data
# df_allelefreq <- fread("data/df_allelefreq_combined.csv")
# df_allelefreq <- df_allelefreq[population != "all"]
# df_allelefreq$frequency <- ifelse(df_allelefreq$frequency == 0, FALLBACK_FREQ, df_allelefreq$frequency)
# df_allelefreq[, allele := as.character(allele)]
# 
# # Load kinship coefficients
# kinship_matrix <- fread("data/kinship_coefficients.csv")
# 
# # Load core loci definitions
# core_loci <- fread("data/core_CODIS_loci.csv")
# loci_lists <- list(
#   core_13 = core_loci[core_13 == 1]$locus,
#   identifiler_15 = core_loci[identifiler_15 == 1]$locus,
#   expanded_20 = core_loci[expanded_20 == 1]$locus
# )
# 
# # Example: Calculate LRs for simulated pairs
# # Assuming pair_data comes from Module 4 (simulate_multiple_pairs)
# tested_relationships <- c("parent_child", "full_siblings", "unrelated")
# tested_populations <- c("AfAm", "Cauc", "Hispanic", "Asian")
# 
# # Calculate single-locus LRs
# single_locus_lr <- calculate_single_locus_lr(
#   pair_data = pair_data,
#   tested_relationship = "full_siblings",
#   tested_populations = tested_populations,
#   allele_frequency_data = df_allelefreq,
#   kinship_coefficients = kinship_matrix
# )
# 
# # Example: Calculate single-locus LRs for one tested relationship
# single_locus_lr <- calculate_single_locus_lr(
#   pair_data = pair_data,
#   tested_relationship = "full_siblings",
#   tested_populations = tested_populations,
#   allele_frequency_data = df_allelefreq,
#   kinship_coefficients = kinship_matrix
# )