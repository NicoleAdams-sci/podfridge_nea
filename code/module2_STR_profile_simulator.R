# ------------------------------------------------------------------------------
# Module 2: STR Profile Simulator
# ------------------------------------------------------------------------------

library(dplyr)
library(data.table)

# Source the utility functions to get FALLBACK_FREQ and other shared constants
source("code/LR_kinship_utility_functions.R")

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
# # Load Allele Frequencies Data
# df_allelefreq <- fread("data/df_allelefreq_combined.csv")
# #df_allelefreq <- df_allelefreq[population != "all"] # Filter out "all" population
# df_allelefreq$frequency <- ifelse(df_allelefreq$frequency == 0, FALLBACK_FREQ, df_allelefreq$frequency) # Use shared constant
# df_allelefreq[, allele := as.character(allele)]
# 
# # Extract unique loci
# loci_list <- df_allelefreq |> pull(marker) |> unique()
# 
# # Load Core Loci Data
#   core_loci <- fread("data/core_CODIS_loci.csv")
#   columns <- c("core_13", "identifiler_15", "expanded_20", "supplementary")
#   loci_lists <- lapply(columns, function(col) {
#     core_loci |>
#       filter(get(col) == 1) |>
#       pull(locus)
#   })
#   names(loci_lists) <- columns
#   loci_lists$autosomal_29 <- loci_list
#   
# # This could be core_13, identifiler_15, or any custom list
# profile <- simulate_str_profile(loci_lists$core_13, "AfAm", df_allelefreq)
