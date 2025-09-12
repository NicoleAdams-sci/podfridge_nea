# ==============================================================================
# Module 1: Allele Simulator
# ==============================================================================
#
# Dependencies:
#   - LR_kinship_utility_functions.R (direct)
#
# Direct function calls:
#   - FALLBACK_FREQ constant from utility functions
#
# ------------------------------------------------------------------------------

library(dplyr)
library(data.table)

# Source the utility functions to get FALLBACK_FREQ and other shared constants
source("code/LR_kinship_utility_functions.R")

# ------------------------------------------------------------------------------
# Module 1: Allele Simulator
# ------------------------------------------------------------------------------
#' Simulate a single allele based on locus, population, and allele frequencies
#' Usage: simulate_allele(locus, population, allele_frequency_data)
#' @param locus Character, the locus name (e.g., "D3S1358")
#' @param population Character, the population (e.g., "AfAm", "Cauc", "Hispanic", "Asian")
#' @param allele_frequency_data Data frame with columns: population, marker, allele, frequency
#' @return Character, a single allele

# Function
# ------------------------------------------------------------------------------
# Module 1: Allele Simulator
# ------------------------------------------------------------------------------
#' Simulate a single allele based on locus, population, and allele frequencies
#'
#' @param locus Character, the locus name (e.g., "D3S1358")
#' @param population Character, the population (e.g., "AfAm", "Cauc", "Hispanic", "Asian")
#' @param allele_frequency_data Data frame with columns: population, marker, allele, frequency
#' @return Character, a single allele

simulate_allele <- function(locus, population, allele_frequency_data) {
  # Filter allele frequencies for the specific population and locus
  # Note: frequency > 0 filter removed since we handle zero frequencies with FALLBACK_FREQ
  allele_freqs <- allele_frequency_data |>
    filter(population == !!population, marker == !!locus)
  
  # Check if valid alleles exist
  if (nrow(allele_freqs) == 0) {
    stop(paste("No valid alleles found for population", population, "and locus", locus))
  }
  
  # Extract alleles and frequencies
  alleles <- allele_freqs$allele
  frequencies <- allele_freqs$frequency
  
  # Apply FALLBACK_FREQ to zero frequencies (consistent with utility functions)
  frequencies <- ifelse(frequencies == 0, FALLBACK_FREQ, frequencies)
  
  # Normalize frequencies to ensure they sum to 1
  frequencies <- frequencies / sum(frequencies)
  
  # Sample one allele based on frequencies
  sampled_allele <- sample(alleles, size = 1, replace = TRUE, prob = frequencies)
  
  return(as.character(sampled_allele))
}


##### Usage
# # Use Module 1 - simulate single allele
# single_allele <- simulate_allele("D3S1358", "AfAm", df_allelefreq)
