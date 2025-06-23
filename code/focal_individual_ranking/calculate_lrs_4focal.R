# Load Required Libraries
suppressMessages(suppressWarnings({
  library(tidyverse)
  library(data.table)
  library(parallel)
  library(future)
  library(furrr)
}))

# Command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript calculate_lrs_4focal.R <population>")
}

population <- args[1]
num_cores <- as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK', unset = 4))

# Helper function for logging
log_message <- function(message) {
  cat(paste0("[", Sys.time(), "] ", message, "\n"))
}

# Helper function to log function timings
timing_log <- list()

log_function_time <- function(func, name, ...) {
  start_time <- Sys.time()
  result <- func(...)
  end_time <- Sys.time()
  duration <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  if (!name %in% names(timing_log)) {
    timing_log[[name]] <- list(total = 0, count = 0, min = Inf, max = -Inf, times = c())
  }
  
  timing_log[[name]]$total <- timing_log[[name]]$total + duration
  timing_log[[name]]$count <- timing_log[[name]]$count + 1
  timing_log[[name]]$min <- min(timing_log[[name]]$min, duration)
  timing_log[[name]]$max <- max(timing_log[[name]]$max, duration)
  timing_log[[name]]$times <- c(timing_log[[name]]$times, duration)
  
  return(result)
}

log_message(paste("Starting LR calculation for population:", population))

# Set up parallel processing
log_message(paste("Setting up parallel processing with", num_cores, "cores"))
cl <- makeCluster(num_cores)
plan(cluster, workers = cl)
on.exit(stopCluster(cl))

# Define paths with correct directory structure from simulation script
focal_db_dir <- file.path("output", "focal_database", population)
combined_relationships_file <- file.path(focal_db_dir, paste0("focal_all_relationships_", population, ".csv"))
allele_freq_file <- "data/df_allelefreq_combined.csv"
output_file <- file.path(focal_db_dir, paste0("focal_all_relationships_with_LR_", population, ".csv"))
timing_log_file <- file.path(focal_db_dir, paste0("timing_log_LR_", population, ".csv"))

# Check if files exist
if (!file.exists(combined_relationships_file)) {
  stop(paste("Error: Combined relationships file not found:", combined_relationships_file))
}

if (!file.exists(allele_freq_file)) {
  stop(paste("Error: Allele frequency file not found:", allele_freq_file))
}

# Create focal database directory if needed
dir.create(focal_db_dir, recursive = TRUE, showWarnings = FALSE)

# Define populations
populations_list <- c("AfAm", "Cauc", "Hispanic", "Asian")

# Define kinship coefficients for ONLY parent_child and full_siblings
# These are the two relationship types we want to calculate LRs for
lr_relationship_types <- c("parent_child", "full_siblings")

kinship_matrix <- tibble(
  relationship_type = factor(
    lr_relationship_types,
    levels = lr_relationship_types
  ),
  k0 = c(0, 1/4),      # parent_child=0, full_siblings=1/4
  k1 = c(1, 1/2),      # parent_child=1, full_siblings=1/2
  k2 = c(0, 1/4)       # parent_child=0, full_siblings=1/4
)

# Load Core Loci Data
log_message("Loading core loci data...")
core_loci_time <- system.time({
  core_loci <- fread("data/core_CODIS_loci.csv")
  columns <- c("core_13", "identifiler_15", "expanded_20", "supplementary")
  loci_lists <- lapply(columns, function(col) {
    core_loci |>
      filter(get(col) == 1) |>
      pull(locus)
  })
  names(loci_lists) <- columns
})
log_message(paste("Loaded core loci data in", core_loci_time["elapsed"], "seconds."))

# Load the combined relationships file from simulation
log_message(paste("Loading relationships data from:", combined_relationships_file))
relationships_data <- fread(combined_relationships_file)
log_message(paste("Loaded", nrow(relationships_data), "rows from the relationships file"))

# Extract unique loci from the database
loci_list <- relationships_data |>
  pull(locus) |>
  unique()

# Add to loci_lists
loci_lists$autosomal_29 <- loci_list

# Load allele frequencies
log_message("Loading allele frequencies...")
allele_freqs <- fread(allele_freq_file)
allele_freqs <- allele_freqs[population != "all"]
allele_freqs$frequency <- ifelse(allele_freqs$frequency == 0, 5/(2*1036), allele_freqs$frequency)
allele_freqs[, allele := as.character(allele)]
log_message(paste("Loaded allele frequencies for", length(unique(allele_freqs$population)), "populations"))

# Original calculate_likelihood_ratio function (SWGDAM approach)
calculate_likelihood_ratio <- function(shared_alleles, genotype_match = NULL, pA = NULL, pB = NULL, k0, k1, k2) {
  if (shared_alleles == 0) {
    LR <- k0
    return(LR)
  }
  if (shared_alleles == 1) {
    if (genotype_match == "AA-AA") {
      Rxp <- pA
    } else if (genotype_match == "AA-AB" | genotype_match == "AB-AA") {
      Rxp <- 2 * pA
    } else if (genotype_match == "AB-AC") {
      Rxp <- 4 * pA
    } else if (genotype_match == "AB-AB") {
      Rxp <- (4 * (pA * pB)) / (pA + pB)
    } else {
      stop("Invalid genotype match for 1 shared allele.")
    }
    LR <- k0 + (k1 / Rxp)
    return(LR)
  }
  if (shared_alleles == 2) {
    if (genotype_match == "AA-AA") {
      Rxp <- pA
      Rxu <- pA^2
    } else if (genotype_match == "AB-AB") {
      Rxp <- (4 * pA * pB) / (pA + pB)
      Rxu <- 2 * pA * pB
    } else {
      stop("Invalid genotype match for 2 shared alleles.")
    }
    LR <- k0 + (k1 / Rxp) + (k2 / Rxu)
    return(LR)
  }
}

# Updated kinship_calculation function for focal data format
kinship_calculation_focal <- function(row, allele_frequency_data, kinship_matrix, population_list, lr_relationship_types) {
  # Use focal_* and ind2_* columns instead of ind1_* and ind2_*
  alleles_focal <- as.character(c(row$focal_allele1, row$focal_allele2))
  alleles_ind2 <- as.character(c(row$ind2_allele1, row$ind2_allele2))
  
  shared_alleles_vector <- intersect(alleles_focal, alleles_ind2)
  unique_alleles_focal <- setdiff(alleles_focal, shared_alleles_vector)
  unique_alleles_ind2 <- setdiff(alleles_ind2, shared_alleles_vector)
  
  allele_map <- list()
  next_label <- 1
  
  for (allele in shared_alleles_vector) {
    allele_map[[LETTERS[next_label]]] <- allele
    next_label <- next_label + 1
  }
  
  for (allele in unique_alleles_focal) {
    if (!(allele %in% allele_map)) {
      allele_map[[LETTERS[next_label]]] <- allele
      next_label <- next_label + 1
    }
  }
  
  for (allele in unique_alleles_ind2) {
    if (!(allele %in% allele_map)) {
      allele_map[[LETTERS[next_label]]] <- allele
      next_label <- next_label + 1
    }
  }
  
  allele_map <- unlist(allele_map)
  labeled_alleles_focal <- sapply(as.character(alleles_focal), function(x) names(allele_map)[which(allele_map == x)])
  labeled_alleles_ind2 <- sapply(as.character(alleles_ind2), function(x) names(allele_map)[which(allele_map == x)])
  
  shared_alleles <- length(shared_alleles_vector)
  genotype_focal <- paste(sort(labeled_alleles_focal), collapse = "")
  genotype_ind2 <- paste(sort(labeled_alleles_ind2), collapse = "")
  genotype_match <- paste(genotype_focal, genotype_ind2, sep = "-")
  
  row$shared_alleles <- shared_alleles
  row$genotype_match <- genotype_match
  
  # Calculate LR for ONLY parent_child and full_siblings relationships, for ALL populations
  A_allele <- ifelse("A" %in% names(allele_map), allele_map[["A"]], NA)
  B_allele <- ifelse("B" %in% names(allele_map), allele_map[["B"]], NA)
  
  # Calculate LR for each relationship type and each population
  for (lr_rel_type in lr_relationship_types) {
    # Get kinship coefficients for this LR relationship type
    k_values <- kinship_matrix[kinship_matrix$relationship_type == lr_rel_type, ]
    
    for (pop in population_list) {
      pop_freqs <- dplyr::filter(allele_frequency_data, population == pop, marker == row$locus)
      
      if (nrow(pop_freqs) > 0) {
        pA <- ifelse(any(pop_freqs$allele == A_allele), pop_freqs$frequency[pop_freqs$allele == A_allele], NA)
        pB <- ifelse(any(pop_freqs$allele == B_allele), pop_freqs$frequency[pop_freqs$allele == B_allele], NA)
        
        if (!is.na(pA)) {
          if (is.na(pB) && length(shared_alleles_vector) > 1) {
            # If we need pB but it's missing, set LR to NA
            row[[paste0("LR_", lr_rel_type, "_", pop)]] <- NA
          } else {
            row[[paste0("LR_", lr_rel_type, "_", pop)]] <- calculate_likelihood_ratio(shared_alleles, genotype_match, pA, pB, k_values$k0, k_values$k1, k_values$k2)
          }
        } else {
          row[[paste0("LR_", lr_rel_type, "_", pop)]] <- NA
        }
      } else {
        row[[paste0("LR_", lr_rel_type, "_", pop)]] <- NA
      }
    }
  }
  
  return(row)
}

# Function to process the database in chunks for memory efficiency
process_focal_database_in_chunks <- function(database, allele_frequency_data, kinship_matrix, population_list, lr_relationship_types, chunk_size = 10000) {
  log_message(paste("Processing focal database in chunks of size", chunk_size))
  
  # Process in chunks
  num_rows <- nrow(database)
  num_chunks <- ceiling(num_rows / chunk_size)
  
  result_database <- data.table()
  
  for (i in 1:num_chunks) {
    start_idx <- (i - 1) * chunk_size + 1
    end_idx <- min(i * chunk_size, num_rows)
    
    log_message(paste("Processing chunk", i, "of", num_chunks, "(rows", start_idx, "to", end_idx, ")"))
    
    # Extract chunk
    chunk <- database[start_idx:end_idx, ]
    
    # Process chunk in parallel
    chunk_results <- chunk %>%
      future_pmap(function(...) {
        row_data <- list(...)
        kinship_calculation_focal(row_data, allele_frequency_data, kinship_matrix, population_list, lr_relationship_types)
      }, .options = furrr_options(seed = TRUE))
    
    # Combine results
    chunk_results_dt <- rbindlist(chunk_results)
    result_database <- rbindlist(list(result_database, chunk_results_dt), fill = TRUE)
    
    # Report progress
    log_message(paste("Completed chunk", i, "of", num_chunks, "(", round(i/num_chunks*100), "% complete)"))
  }
  
  return(result_database)
}

# Calculate LRs for all pairs
log_message("Calculating LRs for all focal relationship pairs...")
lr_calculation_time <- system.time({
  result_database <- log_function_time(
    process_focal_database_in_chunks, 
    "process_focal_database_in_chunks", 
    relationships_data, 
    allele_freqs, 
    kinship_matrix, 
    populations_list,
    lr_relationship_types,
    chunk_size = 10000
  )
})
log_message(paste("LR calculation completed in", lr_calculation_time["elapsed"], "seconds"))

# Save results
log_message(paste("Saving results to:", output_file))
fwrite(result_database, output_file)

# Define function to calculate combined LRs across loci sets for focal data
calculate_combined_lrs_focal <- function(final_results, loci_lists, population_list, lr_relationship_types) {
  final_results <- as.data.table(final_results)
  
  # Create an empty data.table to store our results
  result_dt <- data.table(
    family_id = character(),
    population = character(),
    relationship_type = character(),  # This is the simulated relationship (e.g., "full_siblings_focal")
    pair_id = integer(),
    focal_id = integer(),
    stringsAsFactors = FALSE
  )
  
  # Get unique combinations of family_id, population, relationship_type, pair_id, focal_id
  unique_combinations <- unique(final_results[, .(family_id, population, relationship_type, pair_id, focal_id)])
  
  # Calculate combined LRs for each combination
  for (i in 1:nrow(unique_combinations)) {
    family <- unique_combinations$family_id[i]
    pop <- unique_combinations$population[i]
    rel <- unique_combinations$relationship_type[i]
    pair <- unique_combinations$pair_id[i]
    focal <- unique_combinations$focal_id[i]
    
    # Create a row for this combination
    row_data <- list(
      family_id = family,
      population = pop,
      relationship_type = rel,
      pair_id = pair,
      focal_id = focal
    )
    
    # Filter data for this combination
    subset_data <- final_results[family_id == family & population == pop & relationship_type == rel & pair_id == pair & focal_id == focal]
    
    # Calculate combined LRs for each LR relationship type and population and loci set
    for (lr_rel_type in lr_relationship_types) {
      for (test_pop in population_list) {
        col_name <- paste0("LR_", lr_rel_type, "_", test_pop)
        if (col_name %in% names(subset_data)) {
          for (loci_set_name in names(loci_lists)) {
            loci_set <- loci_lists[[loci_set_name]]
            # Calculate product of LRs for this loci set and population
            lr_values <- subset_data[locus %in% loci_set, get(col_name)]
            lr_values[is.na(lr_values)] <- 1.0  # Neutral LR for missing values
            lr_product <- prod(lr_values, na.rm = TRUE)
            # Add to row data with naming scheme: loci_set_LR_relationship_population
            row_data[[paste0(loci_set_name, "_LR_", lr_rel_type, "_", test_pop)]] <- lr_product
          }
        }
      }
    }
    
    # Convert to data.table and bind to result_dt
    row_dt <- as.data.table(row_data)
    result_dt <- rbindlist(list(result_dt, row_dt), fill = TRUE)
  }
  
  return(result_dt)
}

# Calculate combined LRs
log_message("Calculating combined LRs across loci sets...")
combined_lr_time <- system.time({
  combined_lrs <- log_function_time(
    calculate_combined_lrs_focal,
    "calculate_combined_lrs_focal",
    result_database,
    loci_lists,
    populations_list,
    lr_relationship_types
  )
})
log_message(paste("Combined LR calculation completed in", combined_lr_time["elapsed"], "seconds"))

# Save combined LRs in focal database directory
combined_lr_file <- file.path(focal_db_dir, paste0("focal_combined_LR_", population, ".csv"))
log_message(paste("Saving combined LRs to:", combined_lr_file))
fwrite(combined_lrs, combined_lr_file)

# Create summary statistics for LRs by relationship type
log_message("Creating summary statistics...")
summary_stats <- result_database[, .(
  count = .N,
  families = length(unique(family_id)),
  pairs = length(unique(paste(family_id, pair_id))),
  focal_individuals = length(unique(focal_id))
), by = relationship_type]

log_message("Summary of processed data:")
print(summary_stats)

# Save summary
summary_file <- file.path(focal_db_dir, paste0("focal_processing_summary_", population, ".csv"))
fwrite(summary_stats, summary_file)

# Save timing log to CSV
timing_log_df <- tibble(
  function_name = names(timing_log),
  total_time = sapply(timing_log, function(x) x$total),
  count = sapply(timing_log, function(x) x$count),
  min_time = sapply(timing_log, function(x) x$min),
  max_time = sapply(timing_log, function(x) x$max),
  avg_time = sapply(timing_log, function(x) x$total / x$count)
)
write_csv(timing_log_df, timing_log_file)

log_message("LR calculation process complete!")
log_message(paste("Results saved to:", output_file))
log_message(paste("Combined LRs saved to:", combined_lr_file))
log_message(paste("Summary saved to:", summary_file))