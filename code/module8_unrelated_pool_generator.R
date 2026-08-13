# ------------------------------------------------------------------------------
# Module 8: Unrelated Pool Generator
# ------------------------------------------------------------------------------
#
# Dependencies:
#   - Module 1: Allele Simulator (via Module 2)
#   - Module 2: STR Profile Simulator (direct)
#   - LR_kinship_utility_functions.R (via Module 2)
#
# Direct function calls:
#   - simulate_str_profile() from Module 2
#
# ------------------------------------------------------------------------------
library(dplyr)
library(data.table)
library(purrr)

# Source required modules
source("code/LR_kinship_utility_functions.R")
source("code/module1_allele_simulator.R")
source("code/module2_STR_profile_simulator.R")

#' Generate unrelated individuals for a single population
#'
#' Creates a pool of unrelated individuals for a specific population,
#' then saves results to CSV
#'
#' @param population Character, population for simulation (e.g., "AfAm", "Cauc")
#' @param n_unrelated Integer, number of unrelated individuals to generate
#' @param loci_list Character vector, loci to simulate (defaults to all available)
#' @param allele_frequency_data Data frame with allele frequency data
#' @param output_dir Character, directory to save files (defaults to "output")
#' @param custom_datetime Character, custom datetime string (defaults to current time)
#' @return List with 'data' (the unrelated individuals data frame) and 'file_path' (path to saved CSV)
generate_unrelated_pool <- function(population,
                                     n_unrelated,
                                     loci_list = NULL,
                                     allele_frequency_data,
                                     output_dir = "output",
                                     custom_datetime = NULL) {
  
  # Validate inputs
  if (!population %in% unique(allele_frequency_data$population)) {
    stop(paste("Population", population, "not found in allele_frequency_data"))
  }
  
  if (n_unrelated < 1) {
    stop("n_unrelated must be at least 1")
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
  filename <- paste0(population, "_unrelated_", datetime_str, ".csv")
  file_path <- file.path(output_dir, filename)
  
  cat("Generating", n_unrelated, "unrelated individuals for", population, "...\n")
  
  # Generate unrelated individuals
  unrelated_data <- generate_unrelated_individuals(
    n_unrelated = n_unrelated,
    population = population,
    loci_list = loci_list,
    allele_frequency_data = allele_frequency_data,
    batch_id = datetime_str
  )
  
  # Save to CSV
  cat("Saving results to:", file_path, "\n")
  fwrite(unrelated_data, file_path)
  
  # Print summary
  cat("Successfully generated:\n")
  cat("  - Unrelated individuals:", n_unrelated, "\n")
  cat("  - Population:", population, "\n")
  cat("  - Loci:", length(loci_list), "\n")
  cat("  - Total rows:", nrow(unrelated_data), "\n")
  cat("  - File:", filename, "\n")
  
  return(list(
    data = unrelated_data,
    file_path = file_path
  ))
}

#' Generate unrelated individuals
#'
#' Internal function that creates unrelated individuals.
#'
#' Vectorized rewrite: instead of calling simulate_str_profile()/simulate_allele()
#' once per individual per locus (which re-filters and re-normalizes the whole
#' allele-frequency table on every one of the ~n_unrelated * n_loci * 2 draws),
#' this precomputes the allele pool and normalized probabilities ONCE per locus
#' and then draws all n_unrelated alleles for that locus in a single sample() call.
#'
#' This is identical in distribution to the original: each allele is an
#' independent draw from the same per-(population, locus) categorical
#' distribution (same FALLBACK_FREQ handling, same normalization). A1 and A2 are
#' drawn separately so the two alleles at a locus stay independent, exactly as in
#' simulate_str_profile(). Output columns, column order, individual_id format,
#' and row order (individual-major, loci in loci_list order) are unchanged.
#'
#' @param n_unrelated Integer, number of unrelated individuals
#' @param population Character, population for simulation
#' @param loci_list Character vector, loci to simulate
#' @param allele_frequency_data Data frame with allele frequencies
#' @param batch_id Character, batch identifier
#' @return Data frame with all unrelated individuals
generate_unrelated_individuals <- function(n_unrelated,
                                           population,
                                           loci_list,
                                           allele_frequency_data,
                                           batch_id) {
  
  # Fast per-locus lookups; as.data.table() copies, so the caller's object is
  # never modified.
  freq_dt <- as.data.table(allele_frequency_data)
  target_pop <- population
  
  # --- Precompute allele pool + normalized probabilities ONCE per locus --------
  # Mirrors simulate_allele(): substitute FALLBACK_FREQ for any zero frequency,
  # then normalize to sum to 1. Done a single time per locus instead of once per
  # allele draw.
  locus_tables <- lapply(loci_list, function(this_locus) {
    af <- freq_dt[population == target_pop & marker == this_locus]
    
    if (nrow(af) == 0) {
      stop(paste("No valid alleles found for population", target_pop,
                 "and locus", this_locus))
    }
    
    freqs <- af$frequency
    freqs[freqs == 0] <- FALLBACK_FREQ
    freqs <- freqs / sum(freqs)
    
    list(
      alleles = as.character(af$allele),
      prob    = freqs
    )
  })
  names(locus_tables) <- loci_list
  
  # --- Vectorized draw: all n_unrelated individuals for one locus at once ------
  # sample(size = n_unrelated, replace = TRUE) is n_unrelated independent draws
  # from the same distribution, i.e. distributionally identical to calling
  # simulate_allele() once per individual.
  per_locus <- lapply(loci_list, function(this_locus) {
    tab <- locus_tables[[this_locus]]
    data.table(
      ind_idx = seq_len(n_unrelated),
      locus   = this_locus,
      A1 = sample(tab$alleles, size = n_unrelated, replace = TRUE, prob = tab$prob),
      A2 = sample(tab$alleles, size = n_unrelated, replace = TRUE, prob = tab$prob)
    )
  })
  
  all_individuals <- rbindlist(per_locus)
  
  # Restore individual-major row order (all of individual 1's loci, then
  # individual 2's, ...) with loci in loci_list order within each individual, so
  # the file layout matches the original loop-based output. Ordering on the
  # integer ind_idx (not the padded id string) preserves numeric order for
  # n_unrelated >= 1000.
  all_individuals[, locus := factor(locus, levels = loci_list)]
  setorder(all_individuals, ind_idx, locus)
  all_individuals[, locus := as.character(locus)]
  
  # Attach the constant/derived columns and match the original schema exactly.
  all_individuals[, `:=`(
    batch_id              = batch_id,
    individual_id         = sprintf("unrel_%03d", ind_idx),
    relationship_to_focal = "unrelated",
    population            = target_pop
  )]
  all_individuals[, ind_idx := NULL]
  
  setcolorder(
    all_individuals,
    c("batch_id", "individual_id", "relationship_to_focal",
      "population", "locus", "A1", "A2")
  )
  
  as.data.frame(all_individuals)
}

#' Generate unrelated pools for multiple populations
#'
#' Convenience function to generate unrelated pools for multiple populations
#'
#' @param populations Character vector of populations
#' @param n_unrelated_per_pop Integer, number of unrelated individuals per population
#' @param loci_list Character vector, loci to simulate
#' @param allele_frequency_data Data frame with allele frequency data
#' @param output_dir Character, directory to save files
#' @param use_single_datetime Logical, whether to use same datetime for all files
#' @return Data frame with population details and file paths
generate_multiple_pop_unrelated <- function(populations,
                                             n_unrelated_per_pop,
                                             loci_list = NULL,
                                             allele_frequency_data,
                                             output_dir = "output",
                                             use_single_datetime = TRUE) {
  
  # Generate single datetime if requested
  if (use_single_datetime) {
    shared_datetime <- format(Sys.time(), "%Y%m%d_%H%M%S")
  } else {
    shared_datetime <- NULL
  }
  
  # Generate unrelated pools for each population
  results <- map_dfr(populations, function(pop) {
    result <- generate_unrelated_pool(
      population = pop,
      n_unrelated = n_unrelated_per_pop,
      loci_list = loci_list,
      allele_frequency_data = allele_frequency_data,
      output_dir = output_dir,
      custom_datetime = shared_datetime
    )
    
    data.frame(
      population = pop,
      n_unrelated = n_unrelated_per_pop,
      file_path = result$file_path,
      filename = basename(result$file_path),
      datetime_used = ifelse(is.null(shared_datetime), 
                            format(Sys.time(), "%Y%m%d_%H%M%S"), 
                            shared_datetime),
      stringsAsFactors = FALSE
    )
  })
  
  cat("\nGenerated unrelated pools for", length(populations), "populations.\n")
  
  return(results)
}

# ------------------------------------------------------------------------------
# Usage Examples (commented out)
# ------------------------------------------------------------------------------
# # Single population unrelated pool
# result <- generate_unrelated_pool(
#   population = "AfAm",
#   n_unrelated = 1000,
#   loci_list = loci_list,
#   allele_frequency_data = df_allelefreq
# )
# 
# # Access the data and file path
# unrelated_data <- result$data
# file_path <- result$file_path
# 
# # Multiple populations
# results <- generate_multiple_pop_unrelated(
#   populations = c("AfAm", "Cauc", "Hispanic", "Asian"),
#   n_unrelated_per_pop = 1000,
#   loci_list = loci_list,
#   allele_frequency_data = df_allelefreq
# )
