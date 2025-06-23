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
  stop("Usage: Rscript calculate_lrs_4focal_unrelated.R <population>")
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

log_message(paste("Starting LR calculation for focal-unrelated pairs, population:", population))

# Set up parallel processing
log_message(paste("Setting up parallel processing with", num_cores, "cores"))
cl <- makeCluster(num_cores)
plan(cluster, workers = cl)
on.exit(stopCluster(cl))

# Define paths
focal_db_dir <- file.path("output", "focal_database", population)
input_file <- file.path(focal_db_dir, paste0("focal_unrelated_pairs_", population, ".csv"))
allele_freq_file <- "data/df_allelefreq_combined.csv"
output_file <- file.path(focal_db_dir, paste0("focal_unrelated_with_LR_", population, ".csv"))
combined_output_file <- file.path(focal_db_dir, paste0("focal_unrelated_combined_LR_", population, ".csv"))
timing_log_file <- file.path(focal_db_dir, paste0("timing_log_focal_unrelated_LR_", population, ".csv"))

# Check if files exist
if (!file.exists(input_file)) {
  stop(paste("Error: Focal-unrelated pairs file not found:", input_file))
}

if (!file.exists(allele_freq_file)) {
  stop(paste("Error: Allele frequency file not found:", allele_freq_file))
}

# Create focal database directory if needed
dir.create(focal_db_dir, recursive = TRUE, showWarnings = FALSE)

# Define populations
populations_list <- c("AfAm", "Cauc", "Hispanic", "Asian")

# Define kinship coefficients for parent_child and full_siblings (for unrelated, use unrelated coefficients)
kinship_matrix <- tibble(
  relationship_type = factor(
    c("parent_child", "full_siblings", "unrelated_focal"),
    levels = c("parent_child", "full_siblings", "unrelated_focal")
  ),
  k0 = c(0, 1/4, 1),
  k1 = c(1, 1/2, 0),
  k2 = c(0, 1/4, 0)
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

# Load the focal-unrelated pairs
log_message(paste("Loading focal-unrelated pairs from:", input_file))
focal_unrelated_pairs <- fread(input_file)
log_message(paste("Loaded", nrow(focal_unrelated_pairs), "focal-unrelated pair rows"))

# Check the column names to debug
log_message("Available columns:")
log_message(paste(names(focal_unrelated_pairs), collapse = ", "))

# Extract unique loci from the pairs data
loci_list <- focal_unrelated_pairs |>
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

# Source the LR calculation functions (same as in calculate_lrs_4focal.R)
log_message("Setting up LR calculation functions...")

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

# Kinship calculation function (adapted from original script)
kinship_calculation <- function(row, allele_frequency_data, kinship_matrix, population_list) {
  alleles_ind1 <- as.character(c(row$focal_allele1, row$focal_allele2))
  alleles_ind2 <- as.character(c(row$ind2_allele1, row$ind2_allele2))
  
  shared_alleles_vector <- intersect(alleles_ind1, alleles_ind2)
  unique_alleles_ind1 <- setdiff(alleles_ind1, shared_alleles_vector)
  unique_alleles_ind2 <- setdiff(alleles_ind2, shared_alleles_vector)
  
  allele_map <- list()
  next_label <- 1
  
  for (allele in shared_alleles_vector) {
    allele_map[[LETTERS[next_label]]] <- allele
    next_label <- next_label + 1
  }
  
  for (allele in unique_alleles_ind1) {
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
  
  # Fix the labeling to ensure we get character vectors
  labeled_alleles_ind1 <- character(length(alleles_ind1))
  for (i in 1:length(alleles_ind1)) {
    match_idx <- which(allele_map == alleles_ind1[i])
    if (length(match_idx) > 0) {
      labeled_alleles_ind1[i] <- names(allele_map)[match_idx[1]]
    }
  }
  
  labeled_alleles_ind2 <- character(length(alleles_ind2))
  for (i in 1:length(alleles_ind2)) {
    match_idx <- which(allele_map == alleles_ind2[i])
    if (length(match_idx) > 0) {
      labeled_alleles_ind2[i] <- names(allele_map)[match_idx[1]]
    }
  }
  
  shared_alleles <- length(shared_alleles_vector)
  genotype_ind1 <- paste(sort(labeled_alleles_ind1), collapse = "")
  genotype_ind2 <- paste(sort(labeled_alleles_ind2), collapse = "")
  genotype_match <- paste(genotype_ind1, genotype_ind2, sep = "-")
  
  row$shared_alleles <- shared_alleles
  row$genotype_match <- genotype_match
  
  # Calculate LR for parent_child and full_siblings (treating unrelated_focal as unrelated)
  A_allele <- ifelse("A" %in% names(allele_map), allele_map[["A"]], NA)
  B_allele <- ifelse("B" %in% names(allele_map), allele_map[["B"]], NA)
  
  # Calculate LRs for parent_child and full_siblings hypotheses
  for (rel_type in c("parent_child", "full_siblings")) {
    k_values <- kinship_matrix[kinship_matrix$relationship_type == rel_type, ]
    
    for (pop in population_list) {
      pop_freqs <- dplyr::filter(allele_frequency_data, population == pop, marker == row$locus)
      
      if (nrow(pop_freqs) > 0) {
        pA <- ifelse(any(pop_freqs$allele == A_allele), pop_freqs$frequency[pop_freqs$allele == A_allele], NA)
        pB <- ifelse(any(pop_freqs$allele == B_allele), pop_freqs$frequency[pop_freqs$allele == B_allele], NA)
        
        if (!is.na(pA)) {
          if (is.na(pB) && length(shared_alleles_vector) > 1) {
            # If we need pB but it's missing, set LR to NA
            row[[paste0("LR_", rel_type, "_", pop)]] <- NA
          } else {
            row[[paste0("LR_", rel_type, "_", pop)]] <- calculate_likelihood_ratio(shared_alleles, genotype_match, pA, pB, k_values$k0, k_values$k1, k_values$k2)
          }
        } else {
          row[[paste0("LR_", rel_type, "_", pop)]] <- NA
        }
      } else {
        row[[paste0("LR_", rel_type, "_", pop)]] <- NA
      }
    }
  }
  
  return(row)
}

# Function to process the pairs in chunks for memory efficiency
process_pairs_in_chunks <- function(pairs_data, allele_frequency_data, kinship_matrix, population_list, chunk_size = 10000) {
  log_message(paste("Processing pairs in chunks of size", chunk_size))
  
  # Process in chunks
  num_rows <- nrow(pairs_data)
  num_chunks <- ceiling(num_rows / chunk_size)
  
  result_pairs <- data.table()
  
  for (i in 1:num_chunks) {
    start_idx <- (i - 1) * chunk_size + 1
    end_idx <- min(i * chunk_size, num_rows)
    
    log_message(paste("Processing chunk", i, "of", num_chunks, "(rows", start_idx, "to", end_idx, ")"))
    
    # Extract chunk
    chunk <- pairs_data[start_idx:end_idx, ]
    
    # Process chunk in parallel
    chunk_results <- chunk %>%
      future_pmap(function(...) {
        row_data <- list(...)
        kinship_calculation(row_data, allele_frequency_data, kinship_matrix, population_list)
      }, .options = furrr_options(seed = TRUE))
    
    # Combine results
    chunk_results_dt <- rbindlist(chunk_results)
    result_pairs <- rbindlist(list(result_pairs, chunk_results_dt), fill = TRUE)
    
    # Report progress
    log_message(paste("Completed chunk", i, "of", num_chunks, "(", round(i/num_chunks*100), "% complete)"))
  }
  
  return(result_pairs)
}

# Calculate LRs for all pairs
log_message("Calculating LRs for all focal-unrelated pairs...")
lr_calculation_time <- system.time({
  result_pairs <- log_function_time(
    process_pairs_in_chunks, 
    "process_pairs_in_chunks", 
    focal_unrelated_pairs, 
    allele_freqs, 
    kinship_matrix, 
    populations_list,
    chunk_size = 10000
  )
})
log_message(paste("LR calculation completed in", lr_calculation_time["elapsed"], "seconds"))

# Save results (similar to focal_all_relationships_with_LR_AfAm.csv)
log_message(paste("Saving locus-level results to:", output_file))
fwrite(result_pairs, output_file)

# Calculate combined LRs across loci sets (similar to focal_combined_LR_AfAm.csv)
calculate_combined_lrs <- function(final_results, loci_lists, population_list) {
  final_results <- as.data.table(final_results)
  
  # Create an empty data.table to store our results
  result_dt <- data.table(
    pair_id = integer(),
    stringsAsFactors = FALSE
  )
  
  # Get unique combinations of pair_id (since this is what identifies each focal-unrelated pair)
  unique_combinations <- unique(final_results[, .(pair_id)])
  
  # Calculate combined LRs for each combination
  for (i in 1:nrow(unique_combinations)) {
    pair_id <- unique_combinations$pair_id[i]
    
    # Create a row for this combination
    row_data <- list(
      pair_id = pair_id
    )
    
    # Filter data for this combination
    subset_data <- final_results[pair_id == !!pair_id]
    
    # Calculate combined LRs for each relationship type, loci set, and population
    for (rel_type in c("parent_child", "full_siblings")) {
      for (loci_set_name in names(loci_lists)) {
        loci_set <- loci_lists[[loci_set_name]]
        
        for (pop in population_list) {
          col_name <- paste0("LR_", rel_type, "_", pop)
          
          if (col_name %in% names(subset_data)) {
            # Get LR values for this loci set
            lr_values <- subset_data[locus %in% loci_set, get(col_name)]
            
            # Debug: print first few calculations for troubleshooting
            if (i <= 3 && rel_type == "parent_child" && pop == "AfAm" && loci_set_name == "core_13") {
              log_message(paste("DEBUG: pair_id", pair_id, "loci_set", loci_set_name, "rel_type", rel_type, "pop", pop))
              log_message(paste("LR values:", paste(lr_values, collapse = ", ")))
              log_message(paste("Product:", prod(lr_values, na.rm = TRUE)))
            }
            
            # Keep original values, including 0s and NAs
            lr_product <- prod(lr_values, na.rm = TRUE)
            
            # Add to row data with naming convention matching focal-related pairs
            combined_col_name <- paste0(loci_set_name, "_LR_", rel_type, "_", pop)
            row_data[[combined_col_name]] <- lr_product
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
    calculate_combined_lrs,
    "calculate_combined_lrs",
    result_pairs,
    loci_lists,
    populations_list
  )
})
log_message(paste("Combined LR calculation completed in", combined_lr_time["elapsed"], "seconds"))

# Save combined LRs (similar to focal_combined_LR_AfAm.csv)
log_message(paste("Saving combined LRs to:", combined_output_file))
fwrite(combined_lrs, combined_output_file)

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

# Print summary
log_message("=== SUMMARY ===")
log_message(paste("Population:", population))
log_message(paste("Total focal-unrelated pairs processed:", nrow(result_pairs)))
log_message(paste("Unique focal-unrelated pair combinations:", nrow(combined_lrs)))
log_message(paste("Loci sets analyzed:", paste(names(loci_lists), collapse = ", ")))
log_message(paste("Populations analyzed:", paste(populations_list, collapse = ", ")))
log_message(paste("Relationship types analyzed: parent_child, full_siblings"))

log_message("Focal-unrelated LR calculation process complete!")