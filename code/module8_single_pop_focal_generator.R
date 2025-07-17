# ------------------------------------------------------------------------------
# Module 8: Single Pop Focal Generator
# ------------------------------------------------------------------------------

library(dplyr)
library(data.table)
library(purrr)

# Source required modules
source("code/LR_kinship_utility_functions.R")
source("code/module1_allele_simulator.R")
source("code/module2_STR_profile_simulator.R")
source("code/module4_related_individual_simulator.R")

#' Generate focal individuals with family structures for a single population
#'
#' Creates focal individuals and their related family members according to
#' specified family structure, then saves results to CSV
#'
#' @param population Character, population for simulation (e.g., "AfAm", "Cauc")
#' @param n_focal Integer, number of focal individuals to generate
#' @param family_structure Character, type of family structure (e.g., "nuclear", "siblings_only", "parent_child")
#' @param loci_list Character vector, loci to simulate (defaults to all available)
#' @param allele_frequency_data Data frame with allele frequency data
#' @param kinship_coefficients Data frame with kinship coefficients
#' @param output_dir Character, directory to save files (defaults to "output")
#' @param custom_datetime Character, custom datetime string (defaults to current time)
#' @return List with 'data' (the family data frame) and 'file_path' (path to saved CSV)
generate_single_pop_focal <- function(population,
                                       n_focal,
                                       family_structure,
                                       loci_list = NULL,
                                       allele_frequency_data,
                                       kinship_coefficients,
                                       output_dir = "output",
                                       custom_datetime = NULL) {
  
  # Validate inputs
  if (!population %in% unique(allele_frequency_data$population)) {
    stop(paste("Population", population, "not found in allele_frequency_data"))
  }
  
  if (n_focal < 1) {
    stop("n_focal must be at least 1")
  }
  
  # Define available family structures
  valid_structures <- c("nuclear", "siblings_only", "parent_child", "parent_only", "child_only")
  if (!family_structure %in% valid_structures) {
    stop(paste("family_structure must be one of:", paste(valid_structures, collapse = ", ")))
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
  filename <- paste0(population, "_focal_", datetime_str, ".csv")
  file_path <- file.path(output_dir, filename)
  
  cat("Generating", n_focal, "focal individuals with", family_structure, "structure for", population, "...\n")
  
  # Generate families based on structure
  family_data <- generate_families_by_structure(
    n_focal = n_focal,
    population = population,
    family_structure = family_structure,
    loci_list = loci_list,
    allele_frequency_data = allele_frequency_data,
    kinship_coefficients = kinship_coefficients,
    batch_id = datetime_str
  )
  
  # Save to CSV
  cat("Saving results to:", file_path, "\n")
  fwrite(family_data, file_path)
  
  # Print summary
  n_individuals <- length(unique(family_data$individual_id))
  cat("Successfully generated:\n")
  cat("  - Focal individuals:", n_focal, "\n")
  cat("  - Total individuals:", n_individuals, "\n")
  cat("  - Family structure:", family_structure, "\n")
  cat("  - Loci:", length(loci_list), "\n")
  cat("  - Total rows:", nrow(family_data), "\n")
  cat("  - File:", filename, "\n")
  
  return(list(
    data = family_data,
    file_path = file_path
  ))
}

#' Generate families based on specified structure
#'
#' Internal function that creates family members according to structure type
#'
#' @param n_focal Integer, number of focal individuals
#' @param population Character, population for simulation
#' @param family_structure Character, family structure type
#' @param loci_list Character vector, loci to simulate
#' @param allele_frequency_data Data frame with allele frequencies
#' @param kinship_coefficients Data frame with kinship coefficients
#' @param batch_id Character, batch identifier
#' @return Data frame with all family members
generate_families_by_structure <- function(n_focal,
                                           population,
                                           family_structure,
                                           loci_list,
                                           allele_frequency_data,
                                           kinship_coefficients,
                                           batch_id) {
  
  all_families <- map_dfr(1:n_focal, function(focal_idx) {
    
    # Format focal ID and family ID with leading zeros
    focal_id <- sprintf("%03d", focal_idx)
    family_id <- paste0("fam_", sprintf("%03d", focal_idx))
    
    # Generate focal individual using Module 2
    focal_profile <- simulate_str_profile(loci_list, population, allele_frequency_data)
    focal_profile$individual_id <- paste0("focal_", focal_id)
    focal_profile$relationship_to_focal <- "self"
    focal_profile$focal_id <- focal_id
    focal_profile$family_id <- family_id
    
    # Initialize family with focal individual
    family_members <- list(focal_profile)
    
    # Add family members based on structure
    if (family_structure == "nuclear") {
      # Add parents and siblings
      parent1 <- simulate_related_individual(focal_profile, "parent_child", allele_frequency_data, paste0("parent1_", focal_id))
      parent2 <- simulate_related_individual(focal_profile, "parent_child", allele_frequency_data, paste0("parent2_", focal_id))
      sibling1 <- simulate_related_individual(focal_profile, "full_siblings", allele_frequency_data, paste0("sibling1_", focal_id))
      sibling2 <- simulate_related_individual(focal_profile, "full_siblings", allele_frequency_data, paste0("sibling2_", focal_id))
      
      parent1$focal_id <- focal_id
      parent2$focal_id <- focal_id
      sibling1$focal_id <- focal_id
      sibling2$focal_id <- focal_id
      
      parent1$family_id <- family_id
      parent2$family_id <- family_id
      sibling1$family_id <- family_id
      sibling2$family_id <- family_id
      
      family_members <- c(family_members, list(parent1, parent2, sibling1, sibling2))
      
    } else if (family_structure == "siblings_only") {
      # Add siblings only
      sibling1 <- simulate_related_individual(focal_profile, "full_siblings", allele_frequency_data, paste0("sibling1_", focal_id))
      sibling2 <- simulate_related_individual(focal_profile, "full_siblings", allele_frequency_data, paste0("sibling2_", focal_id))
      
      sibling1$focal_id <- focal_id
      sibling2$focal_id <- focal_id
      
      sibling1$family_id <- family_id
      sibling2$family_id <- family_id
      
      family_members <- c(family_members, list(sibling1, sibling2))
      
    } else if (family_structure == "parent_child") {
      # Add one parent and one child
      parent1 <- simulate_related_individual(focal_profile, "parent_child", allele_frequency_data, paste0("parent_", focal_id))
      child1 <- simulate_related_individual(focal_profile, "parent_child", allele_frequency_data, paste0("child_", focal_id))
      
      parent1$focal_id <- focal_id
      child1$focal_id <- focal_id
      
      parent1$family_id <- family_id
      child1$family_id <- family_id
      
      family_members <- c(family_members, list(parent1, child1))
      
    } else if (family_structure == "parent_only") {
      # Add parents only
      parent1 <- simulate_related_individual(focal_profile, "parent_child", allele_frequency_data, paste0("parent1_", focal_id))
      parent2 <- simulate_related_individual(focal_profile, "parent_child", allele_frequency_data, paste0("parent2_", focal_id))
      
      parent1$focal_id <- focal_id
      parent2$focal_id <- focal_id
      
      parent1$family_id <- family_id
      parent2$family_id <- family_id
      
      family_members <- c(family_members, list(parent1, parent2))
      
    } else if (family_structure == "child_only") {
      # Add children only
      child1 <- simulate_related_individual(focal_profile, "parent_child", allele_frequency_data, paste0("child1_", focal_id))
      child2 <- simulate_related_individual(focal_profile, "parent_child", allele_frequency_data, paste0("child2_", focal_id))
      
      child1$focal_id <- focal_id
      child2$focal_id <- focal_id
      
      child1$family_id <- family_id
      child2$family_id <- family_id
      
      family_members <- c(family_members, list(child1, child2))
    }
    
    # Combine all family members into one data frame
    family_df <- bind_rows(family_members)
    
    # Add batch_id and reorder columns in specified order
    family_df$batch_id <- batch_id
    family_df <- family_df |>
      select(batch_id, family_id, focal_id, individual_id, relationship_to_focal, population, locus, A1, A2)
    
    return(family_df)
  })
  
  return(all_families)
}

#' Generate focal families for multiple populations
#'
#' Convenience function to generate focal families for multiple populations
#'
#' @param populations Character vector of populations
#' @param n_focal_per_pop Integer, number of focal individuals per population
#' @param family_structure Character, family structure type
#' @param loci_list Character vector, loci to simulate
#' @param allele_frequency_data Data frame with allele frequency data
#' @param kinship_coefficients Data frame with kinship coefficients
#' @param output_dir Character, directory to save files
#' @param use_single_datetime Logical, whether to use same datetime for all files
#' @return Data frame with population details and file paths
generate_multiple_pop_focal <- function(populations,
                                         n_focal_per_pop,
                                         family_structure,
                                         loci_list = NULL,
                                         allele_frequency_data,
                                         kinship_coefficients,
                                         output_dir = "output",
                                         use_single_datetime = TRUE) {
  
  # Generate single datetime if requested
  if (use_single_datetime) {
    shared_datetime <- format(Sys.time(), "%Y%m%d_%H%M%S")
  } else {
    shared_datetime <- NULL
  }
  
  # Generate focal families for each population
  results <- map_dfr(populations, function(pop) {
    result <- generate_single_pop_focal(
      population = pop,
      n_focal = n_focal_per_pop,
      family_structure = family_structure,
      loci_list = loci_list,
      allele_frequency_data = allele_frequency_data,
      kinship_coefficients = kinship_coefficients,
      output_dir = output_dir,
      custom_datetime = shared_datetime
    )
    
    data.frame(
      population = pop,
      n_focal = n_focal_per_pop,
      family_structure = family_structure,
      file_path = result$file_path,
      filename = basename(result$file_path),
      datetime_used = ifelse(is.null(shared_datetime), 
                            format(Sys.time(), "%Y%m%d_%H%M%S"), 
                            shared_datetime),
      stringsAsFactors = FALSE
    )
  })
  
  cat("\nGenerated focal families for", length(populations), "populations.\n")
  
  return(results)
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
# # Single population focal generation
# result <- generate_single_pop_focal(
#   population = "AfAm",
#   n_focal = 100,
#   family_structure = "nuclear",
#   loci_list = loci_list,
#   allele_frequency_data = df_allelefreq,
#   kinship_coefficients = kinship_matrix
# )
# 
# # Access the data and file path
# family_data <- result$data
# file_path <- result$file_path
# 
# # Multiple populations
# results <- generate_multiple_pop_focal(
#   populations = c("AfAm", "Cauc", "Hispanic", "Asian"),
#   n_focal_per_pop = 50,
#   family_structure = "siblings_only",
#   loci_list = loci_list,
#   allele_frequency_data = df_allelefreq,
#   kinship_coefficients = kinship_matrix
# )