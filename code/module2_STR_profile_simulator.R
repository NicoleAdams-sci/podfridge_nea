# ------------------------------------------------------------------------------
# Module 2: STR Profile Simulator
# ------------------------------------------------------------------------------
#
# Dependencies:
#   - Module 1: Allele Simulator (direct)
#   - LR_kinship_utility_functions.R (via Module 1)
#
# Direct function calls:
#   - simulate_allele() from Module 1
#
# ------------------------------------------------------------------------------

library(dplyr)
library(data.table)

# Source the utility functions to get FALLBACK_FREQ and other shared constants
source("code/LR_kinship_utility_functions.R")
source("code/module1_allele_simulator.R")

#' Generate complete STR profiles for individuals
#' Usage:
#' @param loci_list Character vector, list of loci to simulate (e.g., c("D3S1358", "TH01", ...))
#' @param population Character, the population (e.g., "AfAm", "Cauc", "Hispanic", "Asian")
#' @param allele_frequency_data Data frame with columns: population, marker, allele, frequency
#' @return Data frame with columns: locus, allele1, allele2
simulate_str_profile <- function(loci_list, population, allele_frequency_data) {
  
  # Initialize results data frame
  profile <- data.frame(
    population = rep(population, length(loci_list)),  # Track which population was used
    locus = loci_list,
    A1 = character(length(loci_list)),
    A2 = character(length(loci_list)),
    stringsAsFactors = FALSE
  )
  
  # Simulate two alleles for each locus (calling allele simulator twice)
  for (i in seq_along(loci_list)) {
    locus <- loci_list[i]
    
    # Call allele simulator twice for each locus
    profile$A1[i] <- simulate_allele(locus, population, allele_frequency_data)
    profile$A2[i] <- simulate_allele(locus, population, allele_frequency_data)
  }
  
  return(profile)
}


##### Usage
# # This could be core_13, identifiler_15, or any custom list
# profile <- simulate_str_profile(loci_lists$core_13, "AfAm", df_allelefreq)
