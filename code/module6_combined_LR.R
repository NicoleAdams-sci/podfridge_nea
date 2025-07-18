# ------------------------------------------------------------------------------
# Module 6: Combined Likelihood Ratio Calculator
# ------------------------------------------------------------------------------

library(dplyr)
library(tibble)
library(purrr)

#' Calculate combined (multi-locus) likelihood ratios across loci sets
#'
#' Takes single-locus LR results from Module 5 and combines them across 
#' different loci sets by multiplying LRs for loci within each set
#'
#' @param single_locus_results Data frame from Module 5 with single-locus LRs
#' @param loci_sets Named list of loci vectors (e.g., core_13, identifiler_15, expanded_20)
#' @return Data frame with combined LRs for each loci set in long format
calculate_combined_lr <- function(single_locus_results, loci_sets) {
  
  # Validate inputs
  required_cols <- c("batch_id", "pair_id", "population", "known_relationship", 
                     "locus", "tested_relationship", "tested_population", "LR")
  missing_cols <- setdiff(required_cols, names(single_locus_results))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  if (!is.list(loci_sets) || is.null(names(loci_sets))) {
    stop("loci_sets must be a named list")
  }
  
  # Calculate combined LRs for each loci set
  combined_results <- map_dfr(names(loci_sets), function(set_name) {
    loci_in_set <- loci_sets[[set_name]]
    
    # Filter to loci in this set and calculate product
    set_results <- single_locus_results |>
      filter(locus %in% loci_in_set) |>
      group_by(batch_id, pair_id, population, known_relationship, 
               tested_relationship, tested_population) |>
      summarize(
        combined_LR = prod(LR, na.rm = TRUE),
        .groups = 'drop'
      ) |>
      mutate(loci_set = set_name)
    
    return(set_results)
  })
  
  # Reorder columns with batch_id first
  combined_results <- combined_results |>
    select(batch_id, pair_id, population, known_relationship, loci_set,
           tested_relationship, tested_population, combined_LR)
  
  return(combined_results)
}


# Note: Use calculate_combined_lr() with explicitly defined loci sets for better control and transparency

# ------------------------------------------------------------------------------
# Usage Example (commented out)
# ------------------------------------------------------------------------------
# # Load core loci definitions
# core_loci <- read.csv("data/core_CODIS_loci.csv")
# loci_lists <- list(
#   core_13 = core_loci[core_loci$core_13 == 1, ]$locus,
#   identifiler_15 = core_loci[core_loci$identifiler_15 == 1, ]$locus,
#   expanded_20 = core_loci[core_loci$expanded_20 == 1, ]$locus
# )
# 
# # Calculate combined LRs (assuming single_locus_lr comes from Module 5)
# combined_lr_results <- calculate_combined_lr(
#   single_locus_results = single_locus_lr,
#   loci_sets = loci_lists
# )
# 
# # Define loci sets explicitly
# core_loci <- read.csv("data/core_CODIS_loci.csv")
# loci_lists <- list(
#   core_13 = core_loci[core_loci$core_13 == 1, ]$locus,
#   identifiler_15 = core_loci[core_loci$identifiler_15 == 1, ]$locus,
#   expanded_20 = core_loci[core_loci$expanded_20 == 1, ]$locus,
#   supplementary = core_loci[core_loci$supplementary == 1, ]$locus
# )
# 
# # Calculate combined LRs with explicit loci sets
# combined_results <- calculate_combined_lr(single_locus_lr, loci_lists)