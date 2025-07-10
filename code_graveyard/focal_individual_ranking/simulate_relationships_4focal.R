# Load Required Libraries
suppressMessages(suppressWarnings({
  library(tidyverse)
  library(data.table)
  library(parallel)
}))

# Command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript simulate_relationships_4focal.R <population> <num_pairs_per_relationship>")
}

population <- args[1]
num_pairs <- as.numeric(args[2])

# Helper function for logging
log_message <- function(message) {
  cat(paste0("[", Sys.time(), "] ", message, "\n"))
}

log_message(paste("Starting relationship simulation for population:", population))
log_message(paste("Number of pairs per relationship type:", num_pairs))

# Define paths with new directory structure
focal_db_dir <- file.path("output", "focal_database", population)
focal_file <- file.path(focal_db_dir, paste0("focal_individual_", population, ".csv"))
allele_freq_file <- "data/df_allelefreq_combined.csv"

# Check if files exist
if (!file.exists(focal_file)) {
  stop(paste("Error: Focal individual file not found:", focal_file))
}

if (!file.exists(allele_freq_file)) {
  stop(paste("Error: Allele frequency file not found:", allele_freq_file))
}

# Create focal database directory if needed
dir.create(focal_db_dir, recursive = TRUE, showWarnings = FALSE)

# Relationship types to simulate with "_focal" suffix
relationship_types <- c("full_siblings_focal", "half_siblings_focal", "cousins_focal", "second_cousins_focal")

# Define kinship coefficients for each relationship (mapping base relationship types)
kinship_matrix <- tibble(
  relationship_type = factor(
    c("parent_child", "full_siblings", "half_siblings", "cousins", "second_cousins", "unrelated"),
    levels = c("parent_child", "full_siblings", "half_siblings", "cousins", "second_cousins", "unrelated")
  ),
  k0 = c(0, 1/4, 1/2, 7/8, 15/16, 1),
  k1 = c(1, 1/2, 1/2, 1/8, 1/16, 0),
  k2 = c(0, 1/4, 0, 0, 0, 0)
)

# Load focal individual data
log_message("Loading focal individual data...")
focal_data <- fread(focal_file)

# Keep original column names as focal_id, focal_allele1, focal_allele2
# No renaming needed - we'll use these as the focal individual

# Get unique focal IDs to create family IDs
unique_focal_ids <- unique(focal_data$focal_id)
log_message(paste("Found", length(unique_focal_ids), "unique focal individuals"))

# Load allele frequencies
log_message("Loading allele frequencies...")
allele_freqs <- fread(allele_freq_file)
allele_freqs <- allele_freqs[population != "all"]
allele_freqs$frequency <- ifelse(allele_freqs$frequency == 0, 5/(2*1036), allele_freqs$frequency)
allele_freqs[, allele := as.character(allele)]

# Function to simulate genotypes based on relationship
simulate_related_genotype <- function(focal_genotype, relationship_focal, pop_allele_freqs) {
  # Extract base relationship type (remove "_focal" suffix)
  base_relationship <- gsub("_focal$", "", relationship_focal)
  
  # Extract information
  locus <- focal_genotype$locus
  focal_alleles <- c(focal_genotype$focal_allele1, focal_genotype$focal_allele2)
  
  # Filter allele frequencies for this locus and population
  locus_freqs <- pop_allele_freqs[
    population == focal_genotype$population & 
    marker == locus,
    .(allele, frequency)
  ]
  
  # Ensure we have allele frequencies
  if (nrow(locus_freqs) == 0) {
    stop(paste("No allele frequencies for locus", locus, "in population", focal_genotype$population))
  }
  
  # Get kinship coefficients using base relationship type
  k_vals <- kinship_matrix[kinship_matrix$relationship_type == base_relationship, ]
  
  # Determine relationship mechanism based on kinship coefficients
  relationship_choice <- sample(
    c('none', 'one', 'both'), 
    size = 1, 
    prob = c(k_vals$k0, k_vals$k1, k_vals$k2)
  )
  
  # Generate simulated individual's alleles based on relationship
  if (relationship_choice == 'none') {
    # No shared alleles
    sim_alleles <- sample(
      locus_freqs$allele, 
      size = 2, 
      replace = TRUE, 
      prob = locus_freqs$frequency
    )
  } else if (relationship_choice == 'one') {
    # One shared allele
    shared_allele <- sample(focal_alleles, size = 1)
    non_shared_allele <- sample(
      locus_freqs$allele, 
      size = 1, 
      replace = TRUE, 
      prob = locus_freqs$frequency
    )
    
    # Randomize the order
    if (runif(1) > 0.5) {
      sim_alleles <- c(shared_allele, non_shared_allele)
    } else {
      sim_alleles <- c(non_shared_allele, shared_allele)
    }
  } else if (relationship_choice == 'both') {
    # Both alleles shared (identical genotype)
    sim_alleles <- focal_alleles
  }
  
  # Create row with simulated genotype
  # Focal individual keeps focal_* labels, simulated individual gets ind2_* labels
  result <- focal_genotype[, .(
    family_id = paste0("family_", focal_genotype$focal_id),  # Add family_id based on focal_id
    population,
    relationship_type = relationship_focal,  # Keep the "_focal" suffix
    pair_id = focal_genotype$focal_id,  # Use focal_id as pair_id
    locus,
    focal_id = focal_genotype$focal_id,
    focal_allele1 = focal_genotype$focal_allele1,
    focal_allele2 = focal_genotype$focal_allele2,
    ind2_allele1 = sim_alleles[1],
    ind2_allele2 = sim_alleles[2]
  )]
  
  return(result)
}

# Simulate pairs for each relationship type
all_simulated_pairs <- list()

for (rel_type in relationship_types) {
  log_message(paste("Simulating", num_pairs, rel_type, "pairs..."))
  
  # Initialize results for this relationship type
  rel_results <- data.table()
  
  # Process each unique focal individual
  for (focal_idx in 1:length(unique_focal_ids)) {
    focal_id <- unique_focal_ids[focal_idx]
    family_id <- paste0("family_", focal_id)
    
    # Create family-specific directory
    family_dir <- file.path(focal_db_dir, family_id)
    dir.create(family_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Get all loci data for this focal individual
    focal_individual_data <- focal_data[focal_data$focal_id == focal_id]
    
    # Create specified number of pairs for this focal individual
    for (pair_idx in 1:num_pairs) {
      # For each locus of the focal individual
      locus_results <- lapply(1:nrow(focal_individual_data), function(i) {
        focal_row <- focal_individual_data[i]
        sim_result <- simulate_related_genotype(focal_row, rel_type, allele_freqs)
        # Override pair_id to create unique pair identifiers within this family
        sim_result$pair_id <- pair_idx
        return(sim_result)
      })
      
      # Combine results for all loci
      pair_results <- rbindlist(locus_results)
      rel_results <- rbind(rel_results, pair_results)
      
      # Save individual pair file in family directory
      pair_output_file <- file.path(family_dir, paste0(rel_type, "_pair_", pair_idx, ".csv"))
      fwrite(pair_results, pair_output_file)
    }
    
    # Save all pairs for this focal individual and relationship type in family directory
    focal_rel_output_file <- file.path(family_dir, paste0(rel_type, "_all_pairs.csv"))
    focal_rel_data <- rel_results[rel_results$family_id == family_id]
    fwrite(focal_rel_data, focal_rel_output_file)
  }
  
  # Add to overall results
  all_simulated_pairs[[rel_type]] <- rel_results
  
  # Save relationship-specific file in main focal database directory
  output_file <- file.path(focal_db_dir, paste0("focal_", rel_type, "_pairs_", population, ".csv"))
  fwrite(rel_results, output_file)
  log_message(paste("Saved", nrow(rel_results), "rows to", output_file))
}

# Combine all simulated pairs
all_pairs <- rbindlist(all_simulated_pairs)

# Save combined file in focal database directory
combined_output_file <- file.path(focal_db_dir, paste0("focal_all_relationships_", population, ".csv"))
fwrite(all_pairs, combined_output_file)
log_message(paste("Saved combined file with", nrow(all_pairs), "rows to", combined_output_file))

# Calculate shared alleles and genotype match for validation
log_message("Calculating shared alleles and genotype match for validation...")

calculate_shared_info <- function(data) {
  data[, `:=`(
    shared_alleles = {
      focal_alleles <- c(focal_allele1, focal_allele2)
      ind2_alleles <- c(ind2_allele1, ind2_allele2)
      length(intersect(focal_alleles, ind2_alleles))
    },
    genotype_match = {
      # Simplified version just for validation
      focal_genotype <- paste(sort(c(focal_allele1, focal_allele2)), collapse = "")
      ind2_genotype <- paste(sort(c(ind2_allele1, ind2_allele2)), collapse = "")
      paste(focal_genotype, ind2_genotype, sep = "-")
    }
  ), by = seq_len(nrow(data))]
  
  return(data)
}

validation_data <- calculate_shared_info(copy(all_pairs))

# Summarize shared alleles by relationship type and family
shared_summary <- validation_data[, .(
  count = .N,
  zero_shared = sum(shared_alleles == 0) / .N * 100,
  one_shared = sum(shared_alleles == 1) / .N * 100,
  two_shared = sum(shared_alleles == 2) / .N * 100
), by = .(family_id, relationship_type)]

# Also create overall summary by relationship type
overall_summary <- validation_data[, .(
  total_families = length(unique(family_id)),
  total_pairs = length(unique(paste(family_id, pair_id))),
  count = .N,
  zero_shared = sum(shared_alleles == 0) / .N * 100,
  one_shared = sum(shared_alleles == 1) / .N * 100,
  two_shared = sum(shared_alleles == 2) / .N * 100
), by = relationship_type]

# Print validation summaries
log_message("Validation summary of shared alleles by relationship type (overall):")
print(overall_summary)

log_message("Validation summary of shared alleles by family and relationship type (first 10 families):")
print(head(shared_summary, 40))  # Show first 10 families (4 relationship types each)

# Save validation summaries in focal database directory
validation_file <- file.path(focal_db_dir, paste0("focal_validation_overall_", population, ".csv"))
fwrite(overall_summary, validation_file)
log_message(paste("Saved overall validation summary to", validation_file))

validation_by_family_file <- file.path(focal_db_dir, paste0("focal_validation_by_family_", population, ".csv"))
fwrite(shared_summary, validation_by_family_file)
log_message(paste("Saved family-specific validation summary to", validation_by_family_file))

log_message("Simulation complete!")