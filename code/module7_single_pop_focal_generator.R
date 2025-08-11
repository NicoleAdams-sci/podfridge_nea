# ------------------------------------------------------------------------------
# Module 7: Single Pop Focal Generator (Flexible Family Structure)
# ------------------------------------------------------------------------------
#
# Dependencies:
#   - Module 1: Allele Simulator (via Modules 2 & 3)
#   - Module 2: STR Profile Simulator (direct)
#   - Module 3: Related Individual Simulator (direct)
#   - LR_kinship_utility_functions.R (via Modules 2 & 3)
#
# Direct function calls:
#   - simulate_str_profile() from Module 2
#   - simulate_related_individual() from Module 3
#
# ------------------------------------------------------------------------------

library(dplyr)
library(data.table)
library(purrr)

# Source required modules
source("code/LR_kinship_utility_functions.R")
source("code/module1_allele_simulator.R")
source("code/module2_STR_profile_simulator.R")
source("code/module3_related_individual_simulator.R")

#' Generate focal individuals with flexible family structures for a single population
#'
#' Creates focal individuals and their specified relatives according to
#' user-defined family composition, then saves results to CSV
#'
#' @param population Character, population for simulation (e.g., "AfAm", "Cauc")
#' @param n_focal Integer, number of focal individuals to generate
#' @param relationship_counts Named list specifying number of each relationship type.
#'        Valid relationship types: "parent_child", "full_siblings", "half_siblings", 
#'        "first_cousins", "second_cousins", "unrelated"
#'        Example: list(parent_child = 2, full_siblings = 3, first_cousins = 1)
#' @param loci_list Character vector, loci to simulate (defaults to all available)
#' @param allele_frequency_data Data frame with allele frequency data
#' @param kinship_coefficients Data frame with kinship coefficients
#' @param output_dir Character, directory to save files (defaults to "output")
#' @param custom_datetime Character, custom datetime string (defaults to current time)
#' @return List with 'data' (the family data frame) and 'file_path' (path to saved CSV)
generate_single_pop_focal <- function(population,
                                       n_focal,
                                       relationship_counts,
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
  
  # Validate relationship_counts
  if (!is.list(relationship_counts) || is.null(names(relationship_counts))) {
    stop("relationship_counts must be a named list")
  }
  
  # Define valid relationship types (must match kinship_coefficients)
  valid_relationships <- kinship_coefficients$relationship_type
  invalid_relationships <- setdiff(names(relationship_counts), valid_relationships)
  if (length(invalid_relationships) > 0) {
    stop(paste("Invalid relationship types:", paste(invalid_relationships, collapse = ", "),
               "\nValid types:", paste(valid_relationships, collapse = ", ")))
  }
  
  # Validate counts are positive integers
  if (any(unlist(relationship_counts) < 0) || any(!sapply(relationship_counts, is.numeric))) {
    stop("All relationship counts must be non-negative integers")
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
  
  # Generate descriptive filename based on relationship composition
  rel_summary <- paste(names(relationship_counts), unlist(relationship_counts), sep = "", collapse = "_")
  filename <- paste0(population, "_focal_", rel_summary, "_", datetime_str, ".csv")
  file_path <- file.path(output_dir, filename)
  
  cat("Generating", n_focal, "focal individuals for", population, "with relationship structure:\n")
  for (rel in names(relationship_counts)) {
    cat("  -", relationship_counts[[rel]], rel, "\n")
  }
  
  # Generate families based on relationship counts
  family_data <- generate_families_by_relationships(
    n_focal = n_focal,
    population = population,
    relationship_counts = relationship_counts,
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
  total_relatives <- sum(unlist(relationship_counts))
  cat("Successfully generated:\n")
  cat("  - Focal individuals:", n_focal, "\n")
  cat("  - Relatives per focal:", total_relatives, "\n")
  cat("  - Total individuals:", n_individuals, "\n")
  cat("  - Loci:", length(loci_list), "\n")
  cat("  - Total rows:", nrow(family_data), "\n")
  cat("  - File:", basename(filename), "\n")
  
  return(list(
    data = family_data,
    file_path = file_path
  ))
}

#' Generate families based on specified relationship counts
#'
#' Internal function that creates family members according to relationship specifications
#'
#' @param n_focal Integer, number of focal individuals
#' @param population Character, population for simulation
#' @param relationship_counts Named list of relationship counts
#' @param loci_list Character vector, loci to simulate
#' @param allele_frequency_data Data frame with allele frequencies
#' @param kinship_coefficients Data frame with kinship coefficients
#' @param batch_id Character, batch identifier
#' @return Data frame with all family members
generate_families_by_relationships <- function(n_focal,
                                               population,
                                               relationship_counts,
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
    
    # Generate relatives for each relationship type
    for (relationship_type in names(relationship_counts)) {
      n_relatives <- relationship_counts[[relationship_type]]
      
      if (n_relatives > 0) {
        # Generate each relative of this type
        for (rel_idx in 1:n_relatives) {
          
          # Create unique individual ID
          individual_id <- paste0(gsub("_", "", relationship_type), rel_idx, "_", focal_id)
          
          # Generate related individual using Module 3
          related_individual <- simulate_related_individual(
            focal_profile = focal_profile,
            known_relationship = relationship_type,
            allele_frequency_data = allele_frequency_data,
            individual_id = individual_id
          )
          
          # Add family metadata
          related_individual$focal_id <- focal_id
          related_individual$family_id <- family_id
          
          family_members <- c(family_members, list(related_individual))
        }
      }
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

#' Generate focal families for multiple populations with different structures
#'
#' Convenience function to generate focal families for multiple populations,
#' each potentially with different family structures
#'
#' @param population_structures List where each element is named by population
#'        and contains the relationship_counts for that population
#'        Example: list(AfAm = list(parent_child = 2, full_siblings = 4),
#'                     Cauc = list(parent_child = 2, full_siblings = 2))
#' @param n_focal_per_pop Integer, number of focal individuals per population
#' @param loci_list Character vector, loci to simulate
#' @param allele_frequency_data Data frame with allele frequency data
#' @param kinship_coefficients Data frame with kinship coefficients
#' @param output_dir Character, directory to save files
#' @param use_single_datetime Logical, whether to use same datetime for all files
#' @return Data frame with population details and file paths
generate_multiple_pop_focal <- function(population_structures,
                                         n_focal_per_pop,
                                         loci_list = NULL,
                                         allele_frequency_data,
                                         kinship_coefficients,
                                         output_dir = "output",
                                         use_single_datetime = TRUE) {
  
  # Validate input structure
  if (!is.list(population_structures) || is.null(names(population_structures))) {
    stop("population_structures must be a named list")
  }
  
  # Generate single datetime if requested
  if (use_single_datetime) {
    shared_datetime <- format(Sys.time(), "%Y%m%d_%H%M%S")
  } else {
    shared_datetime <- NULL
  }
  
  # Generate focal families for each population
  results <- map_dfr(names(population_structures), function(pop) {
    relationship_counts <- population_structures[[pop]]
    
    result <- generate_single_pop_focal(
      population = pop,
      n_focal = n_focal_per_pop,
      relationship_counts = relationship_counts,
      loci_list = loci_list,
      allele_frequency_data = allele_frequency_data,
      kinship_coefficients = kinship_coefficients,
      output_dir = output_dir,
      custom_datetime = shared_datetime
    )
    
    data.frame(
      population = pop,
      n_focal = n_focal_per_pop,
      relationship_structure = paste(names(relationship_counts), unlist(relationship_counts), sep = "", collapse = "_"),
      file_path = result$file_path,
      filename = basename(result$file_path),
      datetime_used = ifelse(is.null(shared_datetime), 
                            format(Sys.time(), "%Y%m%d_%H%M%S"), 
                            shared_datetime),
      stringsAsFactors = FALSE
    )
  })
  
  cat("\nGenerated focal families for", length(population_structures), "populations.\n")
  
  return(results)
}

#' Create standard relationship structures for common family types
#'
#' Helper function to create relationship_counts for typical family structures
#'
#' @param family_type Character, one of: "nuclear", "large_family", "extended", "minimal"
#' @return Named list with relationship counts
create_standard_family_structure <- function(family_type) {
  
  structures <- list(
    nuclear = list(parent_child = 2, full_siblings = 2),
    large_family = list(parent_child = 2, full_siblings = 4, first_cousins = 2),
    extended = list(parent_child = 2, full_siblings = 3, half_siblings = 1, first_cousins = 3, second_cousins = 2),
    minimal = list(parent_child = 1, full_siblings = 1),
    siblings_only = list(full_siblings = 3),
    parent_only = list(parent_child = 2)
  )
  
  if (!family_type %in% names(structures)) {
    stop(paste("Unknown family_type. Available types:", paste(names(structures), collapse = ", ")))
  }
  
  return(structures[[family_type]])
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
# # Example 1: Custom family structure
# custom_structure <- list(
#   parent_child = 2,      # 2 parents
#   full_siblings = 4,     # 4 full siblings
#   first_cousins = 2      # 2 first cousins
# )
# 
# result1 <- generate_single_pop_focal(
#   population = "AfAm",
#   n_focal = 50,
#   relationship_counts = custom_structure,
#   loci_list = loci_list,
#   allele_frequency_data = df_allelefreq,
#   kinship_coefficients = kinship_matrix
# )
# 
# # Example 2: Different structures for different populations
# population_structures <- list(
#   AfAm = list(parent_child = 2, full_siblings = 4, first_cousins = 3),  # Large families
#   Cauc = list(parent_child = 2, full_siblings = 2),                     # Nuclear families
#   Hispanic = list(parent_child = 2, full_siblings = 3, first_cousins = 2), # Extended families
#   Asian = list(parent_child = 2, full_siblings = 1)                     # Small families
# )
# 
# results2 <- generate_multiple_pop_focal(
#   population_structures = population_structures,
#   n_focal_per_pop = 100,
#   loci_list = loci_list,
#   allele_frequency_data = df_allelefreq,
#   kinship_coefficients = kinship_matrix
# )
# 
# # Example 3: Using predefined family structures
# nuclear_structure <- create_standard_family_structure("nuclear")
# result3 <- generate_single_pop_focal(
#   population = "Cauc",
#   n_focal = 25,
#   relationship_counts = nuclear_structure,
#   loci_list = loci_list,
#   allele_frequency_data = df_allelefreq,
#   kinship_coefficients = kinship_matrix
# )