# Load Required Libraries
suppressMessages(suppressWarnings({
  library(tidyverse)
  library(data.table)
  library(parallel)
}))

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

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript code/sims_post_analysis_NEA.R <input_dir> <population> [output_dir]")
}

input_dir <- args[1]
target_population <- args[2]

# If output directory is provided, use it, otherwise create one based on the population
if (length(args) >= 3) {
  output_dir <- args[3]
} else {
  timestamp <- format(Sys.time(), "%Y%m%d")
  output_dir <- file.path("output", paste0("summary_", target_population, "_", timestamp))
}

# Create output directory if it doesn't exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Define output file paths
combined_lrs_file <- file.path(output_dir, paste0("sim_summary_genotypes_", target_population, ".csv"))
cutoffs_file <- file.path(output_dir, paste0("sim_cutoffs_", target_population, ".csv"))
proportions_file <- file.path(output_dir, paste0("sim_proportions_exceeding_cutoffs_", target_population, ".csv"))
summary_stats_file <- file.path(output_dir, paste0("sim_lr_summary_stats_", target_population, ".csv"))
ratio_summary_file <- file.path(output_dir, paste0("sim_lr_ratio_summary_", target_population, ".csv"))
timing_log_file <- file.path(output_dir, paste0("timing_log_summary_", target_population, ".csv"))

# Define relationship types
relationship_types <- c("parent_child", "full_siblings", "half_siblings", "cousins", "second_cousins", "unrelated")

# Define populations
populations_list <- c("AfAm", "Cauc", "Hispanic", "Asian")

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
  
  # We'll set autosomal_29 when we load the data and get the unique loci
})
log_message(paste("Loaded core loci data in", core_loci_time["elapsed"], "seconds."))

# Load and combine simulation results for the specified population
log_message(paste("Loading simulation results for population:", target_population))
load_time <- system.time({
  final_results <- data.table()
  
  for (rel in relationship_types) {
    file_path <- file.path(input_dir, paste0("sim_processed_genotypes_", target_population, "_", rel, "_combined.csv"))
    
    if (file.exists(file_path)) {
      log_message(paste("Loading file:", file_path))
      data <- fread(file_path)
      final_results <- rbindlist(list(final_results, data), fill = TRUE)
    } else {
      log_message(paste("Warning: File not found:", file_path))
    }
  }
})
log_message(paste("Loaded simulation results in", load_time["elapsed"], "seconds."))

# Check if we have data
if (nrow(final_results) == 0) {
  stop(paste("No data was loaded for population:", target_population))
}

# Extract unique loci from the loaded data
loci_list <- final_results |>
  pull(locus) |>
  unique()

# Add to loci_lists
loci_lists$autosomal_29 <- loci_list

# Define functions for calculations
calculate_combined_lrs <- function(final_results, loci_lists, population_list = populations_list) {
  final_results <- as.data.table(final_results)
  
  # Create an empty data.table to store our results
  result_dt <- data.table(
    population = character(),
    relationship_type = character(),
    sim_id = integer(),
    stringsAsFactors = FALSE
  )
  
  # Get unique combinations of population, relationship_type, and sim_id
  unique_combinations <- unique(final_results[, .(population, relationship_type, sim_id)])
  
  # Calculate combined LRs for each combination
  for (i in 1:nrow(unique_combinations)) {
    pop <- unique_combinations$population[i]
    rel <- unique_combinations$relationship_type[i]
    sim <- unique_combinations$sim_id[i]
    
    # Create a row for this combination
    row_data <- list(
      population = pop,
      relationship_type = rel,
      sim_id = sim
    )
    
    # Filter data for this combination
    subset_data <- final_results[population == pop & relationship_type == rel & sim_id == sim]
    
    # Calculate combined LRs for each population and loci set
    for (test_pop in population_list) {
      col_name <- paste0("LR_", test_pop)
      if (col_name %in% names(subset_data)) {
        for (loci_set_name in names(loci_lists)) {
          loci_set <- loci_lists[[loci_set_name]]
          # Calculate product of LRs for this loci set and population
          lr_product <- prod(subset_data[locus %in% loci_set, get(col_name)], na.rm = TRUE)
          # Add to row data
          row_data[[paste0(loci_set_name, "_", test_pop)]] <- lr_product
        }
      }
    }
    
    # Convert to data.table and bind to result_dt
    row_dt <- as.data.table(row_data)
    result_dt <- rbindlist(list(result_dt, row_dt), fill = TRUE)
  }
  
  # Melt the data
  id_vars <- c("population", "relationship_type", "sim_id")
  measure_vars <- setdiff(names(result_dt), id_vars)
  
  combined_lrs <- melt(result_dt,
                       id.vars = id_vars,
                       measure.vars = measure_vars,
                       variable.name = "loci_set", 
                       value.name = "LR")
  
  # Parse the loci_set to extract base loci set and frequency population
  #combined_lrs[, c("base_loci_set", "freq_population") := tstrsplit(as.character(loci_set), "_", fixed = TRUE, keep = c(1, 2))] # breaks for loci sets with numbers 
  # Add a column to identify if this is the correct population match
  #combined_lrs[, is_correct_pop := (population == freq_population)]
  combined_lrs <- combined_lrs %>%
    mutate(
      freq_population = str_extract(loci_set, "[^_]+$"),
      base_loci_set = str_remove(loci_set, "_[^_]+$")
    ) %>%
    mutate(
      is_correct_pop = (population == freq_population)
    )
  
  
  return(combined_lrs)
}

# Calculate summary statistics for LRs
calculate_summary_stats <- function(combined_lrs) {
  log_message("Calculating summary statistics...")
  summary_stats <- combined_lrs |>
    group_by(relationship_type, population, base_loci_set, freq_population, is_correct_pop) |>
    summarize(
      n = n(),
      mean_LR = mean(LR, na.rm = TRUE),
      median_LR = median(LR, na.rm = TRUE),
      sd_LR = sd(LR, na.rm = TRUE),
      min_LR = min(LR, na.rm = TRUE),
      max_LR = max(LR, na.rm = TRUE),
      lower_95 = quantile(LR, 0.025, na.rm = TRUE),
      upper_95 = quantile(LR, 0.975, na.rm = TRUE),
      .groups = 'drop'
    ) |>
    ungroup()
  
  return(summary_stats)
}

# Calculate ratio of wrong population LR to correct population LR
calculate_ratio_stats <- function(combined_lrs) {
  log_message("Calculating ratio statistics...")
  
  # Extract correct population LRs
  combined_lrs_correct <- combined_lrs[is_correct_pop == TRUE, .(
    population, 
    relationship_type, 
    sim_id, 
    base_loci_set, 
    correct_LR = LR
  )]
  
  # Extract wrong population LRs
  combined_lrs_wrong <- combined_lrs[is_correct_pop == FALSE, .(
    population, 
    relationship_type, 
    sim_id, 
    base_loci_set, 
    freq_population,
    wrong_LR = LR
  )]
  
  # Merge correct and wrong LRs
  combined_lrs_ratio <- merge(
    combined_lrs_wrong,
    combined_lrs_correct,
    by = c("population", "relationship_type", "sim_id", "base_loci_set")
  )
  
  # Calculate ratio of wrong to correct LR
  combined_lrs_ratio[, ratio := wrong_LR / correct_LR]
  
  # Calculate summary statistics for ratios
  ratio_summary <- combined_lrs_ratio |>
    group_by(relationship_type, population, base_loci_set, freq_population) |>
    summarize(
      n = n(),
      mean_ratio = mean(ratio, na.rm = TRUE),
      median_ratio = median(ratio, na.rm = TRUE),
      sd_ratio = sd(ratio, na.rm = TRUE),
      min_ratio = min(ratio, na.rm = TRUE),
      max_ratio = max(ratio, na.rm = TRUE),
      lower_95 = quantile(ratio, 0.025, na.rm = TRUE),
      upper_95 = quantile(ratio, 0.975, na.rm = TRUE),
      .groups = 'drop'
    ) |>
    ungroup()
  
  return(ratio_summary)
}

# Function to calculate cut-off values for 1%, 0.1%, and 0.01% FPR
calculate_cutoffs <- function(input_df, fp_rates) {
  cutoffs <- input_df |>
    filter(relationship_type == "unrelated") |>
    group_by(base_loci_set, freq_population) |>
    summarize(
      fixed_cutoff = 1.00,
      cutoff_1 = quantile(LR, probs = 1 - fp_rates[1] / 100, na.rm = TRUE),
      cutoff_0_1 = quantile(LR, probs = 1 - fp_rates[2] / 100, na.rm = TRUE),
      cutoff_0_01 = quantile(LR, probs = 1 - fp_rates[3] / 100, na.rm = TRUE),
      n_unrelated = n(),
      .groups = 'drop'
    )
  return(cutoffs)
}

calculate_proportions_exceeding_cutoffs <- function(input_df, cutoffs) {
  # Join with cutoffs
  df_with_cutoffs <- left_join(
    input_df,
    cutoffs,
    by = c("base_loci_set", "freq_population")
  )
  
  df_with_cutoffs <- df_with_cutoffs |>
    mutate(
      exceeds_fixed_cutoff = LR > fixed_cutoff,
      exceeds_cutoff_1 = LR > cutoff_1,
      exceeds_cutoff_0_1 = LR > cutoff_0_1,
      exceeds_cutoff_0_01 = LR > cutoff_0_01
    )
  
  proportions_exceeding <- df_with_cutoffs |>
    group_by(population, relationship_type, base_loci_set, freq_population, is_correct_pop) |>
    summarize(
      proportion_exceeding_fixed = sum(exceeds_fixed_cutoff, na.rm = TRUE) / n(),
      proportion_exceeding_1 = sum(exceeds_cutoff_1, na.rm = TRUE) / n(),
      proportion_exceeding_0_1 = sum(exceeds_cutoff_0_1, na.rm = TRUE) / n(),
      proportion_exceeding_0_01 = sum(exceeds_cutoff_0_01, na.rm = TRUE) / n(),
      n_related = n(),
      .groups = 'drop'
    ) |>
    filter(relationship_type != "unrelated")
  
  return(proportions_exceeding)
}

# MAIN EXECUTION SECTION
log_message(paste("Starting summary analysis for population:", target_population))

# Calculate combined LRs
combined_lrs <- log_function_time(calculate_combined_lrs, "calculate_combined_lrs", final_results, loci_lists, populations_list)
fwrite(combined_lrs, combined_lrs_file)
log_message(paste("Combined LRs saved to", combined_lrs_file))

# Calculate and save summary statistics
summary_stats <- log_function_time(calculate_summary_stats, "calculate_summary_stats", combined_lrs)
fwrite(summary_stats, summary_stats_file)
log_message(paste("Summary statistics saved to", summary_stats_file))

# Calculate and save ratio statistics
ratio_summary <- log_function_time(calculate_ratio_stats, "calculate_ratio_stats", combined_lrs)
fwrite(ratio_summary, ratio_summary_file)
log_message(paste("Ratio statistics saved to", ratio_summary_file))

# Calculate and save cutoffs
cutoffs <- log_function_time(calculate_cutoffs, "calculate_cutoffs", combined_lrs, c(1, 0.1, 0.01))
fwrite(cutoffs, cutoffs_file)
log_message(paste("Cutoffs saved to", cutoffs_file))

# Calculate and save proportions exceeding cutoffs
proportions_exceeding_cutoffs <- log_function_time(
  calculate_proportions_exceeding_cutoffs,
  "calculate_proportions_exceeding_cutoffs",
  combined_lrs,
  cutoffs
)
fwrite(proportions_exceeding_cutoffs, proportions_file)
log_message(paste("Proportions exceeding cutoffs saved to", proportions_file))

# Create timing log dataframe
timing_log_df <- tibble(
  function_name = names(timing_log),
  total_time = sapply(timing_log, function(x) x$total),
  count = sapply(timing_log, function(x) x$count),
  min_time = sapply(timing_log, function(x) x$min),
  max_time = sapply(timing_log, function(x) x$max),
  avg_time = sapply(timing_log, function(x) x$total / x$count)
)

# Save timing log
write_csv(timing_log_df, timing_log_file)
log_message(paste("Timing log saved to", timing_log_file))

# Final message
log_message(paste("All summary analyses complete for population:", target_population))