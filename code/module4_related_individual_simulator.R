# ------------------------------------------------------------------------------
# Module 4: Related Individual Simulator
# ------------------------------------------------------------------------------

library(dplyr)
library(data.table)

# Source the utility functions to get FALLBACK_FREQ and other shared constants
source("code/LR_kinship_utility_functions.R")

#' Calculate probability of sharing alleles based on relationship
#'
#' @param relationship_type Character, the known relationship
#' @param kinship_coefficients Data frame with relationship_type, k0, k1, k2
#' @return Character, one of "none", "one", "both" (0, 1, or 2 shared alleles)
calculate_shared_allele_probability <- function(relationship_type, kinship_coefficients) {
  
  # Get kinship coefficients for this relationship
  kinship_coeffs <- kinship_coefficients |>
    filter(relationship_type == !!relationship_type)
  
  if (nrow(kinship_coeffs) != 1) {
    stop(paste("Relationship type", relationship_type, "not found in kinship coefficients"))
  }
  
  # Sample sharing pattern based on k0, k1, k2 probabilities
  sharing_choice <- sample(
    c('none', 'one', 'both'), 
    size = 1, 
    prob = c(kinship_coeffs$k0, kinship_coeffs$k1, kinship_coeffs$k2)
  )
  
  return(sharing_choice)
}


#' Simulate a related individual based on focal profile and known relationship
#'
#' @param focal_profile Data frame with columns: population, locus, A1, A2
#' @param known_relationship Character, the relationship type to simulate
#' @param allele_frequency_data Data frame with population, marker, allele, frequency
#' @param individual_id Character, identifier for the simulated individual
#' @return Data frame with simulated related individual profile
simulate_related_individual <- function(focal_profile, known_relationship, allele_frequency_data, individual_id = "related_ind") {
  
  # Validate inputs
  if (!all(c("population", "locus", "A1", "A2") %in% names(focal_profile))) {
    stop("focal_profile must have columns: population, locus, A1, A2")
  }
  
  if (!known_relationship %in% kinship_matrix$relationship_type) {
    stop(paste("Invalid relationship. Must be one of:", 
               paste(kinship_matrix$relationship_type, collapse = ", ")))
  }
  
  # Initialize result with same structure as focal profile
  related_profile <- focal_profile |>
    select(population, locus) |>
    mutate(
      A1 = character(nrow(focal_profile)),
      A2 = character(nrow(focal_profile)),
      individual_id = individual_id,
      relationship_to_focal = known_relationship
    )
  
  # Simulate alleles for each locus
  for (i in 1:nrow(focal_profile)) {
    focal_row <- focal_profile[i, ]
    
    # Get focal individual's alleles
    focal_alleles <- c(focal_row$A1, focal_row$A2)
    
    # Determine sharing pattern for this locus based on relationship
    sharing_pattern <- calculate_shared_allele_probability(known_relationship, kinship_matrix)
    
    # Simulate related individual's alleles based on sharing pattern
    if (sharing_pattern == 'none') {
      # Share no alleles - simulate independently
      related_profile$A1[i] <- simulate_allele(focal_row$locus, focal_row$population, allele_frequency_data)
      related_profile$A2[i] <- simulate_allele(focal_row$locus, focal_row$population, allele_frequency_data)
      
    } else if (sharing_pattern == 'one') {
      # Share one allele
      shared_allele <- sample(focal_alleles, size = 1)
      non_shared_allele <- simulate_allele(focal_row$locus, focal_row$population, allele_frequency_data)
      
      # Randomly assign which position gets the shared allele
      if (runif(1) > 0.5) {
        related_profile$A1[i] <- shared_allele
        related_profile$A2[i] <- non_shared_allele
      } else {
        related_profile$A1[i] <- non_shared_allele
        related_profile$A2[i] <- shared_allele
      }
      
    } else if (sharing_pattern == 'both') {
      # Share both alleles (same genotype)
      related_profile$A1[i] <- focal_row$A1
      related_profile$A2[i] <- focal_row$A2
    }
  }
  
  # Reorder columns for consistency
  related_profile <- related_profile[, c("individual_id", "population", "locus", "A1", "A2", "relationship_to_focal")]
  
  return(related_profile)
}

# ------------------------------------------------------------------------------
# Wrapper function for pair simulation
# ------------------------------------------------------------------------------
#' Simulate a pair of individuals (focal + related) in wide format
#'
#' @param loci_list Character vector, loci to simulate
#' @param population Character, population for simulation
#' @param known_relationship Character, relationship between individuals
#' @param allele_frequency_data Data frame with allele frequencies
#' @param pair_id Character or numeric, unique identifier for this pair
#' @return Data frame with pair in wide format (one row per locus)
simulate_individual_pair <- function(loci_list, population, known_relationship, allele_frequency_data, pair_id = 1) {
  
  # Step 1: Simulate focal individual de novo (using Module 2)
  focal_profile <- simulate_str_profile(loci_list, population, allele_frequency_data)
  
  # Step 2: Simulate related individual based on focal
  related_profile <- simulate_related_individual(
    focal_profile, 
    known_relationship, 
    allele_frequency_data, 
    "related"
  )
  
  # Step 3: Combine into wide format (one row per locus)
  pair_data <- data.frame(
    pair_id = rep(pair_id, length(loci_list)),
    population = rep(population, length(loci_list)),
    known_relationship = rep(known_relationship, length(loci_list)),
    locus = loci_list,
    focal_A1 = focal_profile$A1,
    focal_A2 = focal_profile$A2,
    ind2_A1 = related_profile$A1[related_profile$individual_id == "related"],
    ind2_A2 = related_profile$A2[related_profile$individual_id == "related"],
    stringsAsFactors = FALSE
  )
  
  return(pair_data)
}

# ------------------------------------------------------------------------------
# Multiple pairs simulation function
# ------------------------------------------------------------------------------
#' Simulate multiple pairs of individuals
#'
#' @param n_pairs Integer, number of pairs to simulate
#' @param loci_list Character vector, loci to simulate
#' @param population Character, population for simulation
#' @param known_relationship Character, relationship between individuals
#' @param allele_frequency_data Data frame with allele frequencies
#' @return Data frame with all pairs in wide format
simulate_multiple_pairs <- function(n_pairs, loci_list, population, known_relationship, allele_frequency_data) {
  
  # Generate all pairs
  pair_list <- lapply(1:n_pairs, function(i) {
    simulate_individual_pair(loci_list, population, known_relationship, allele_frequency_data, pair_id = i)
  })
  
  # Combine all pairs
  all_pairs <- do.call(rbind, pair_list)
  
  return(all_pairs)
}