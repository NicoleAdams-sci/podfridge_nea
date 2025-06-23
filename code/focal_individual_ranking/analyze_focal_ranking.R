# Focal Individual Ranking Analysis (SWGDAM approach) - Individual Pair Testing
library(tidyverse)
library(data.table)

# Command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  population <- "AfAm"  # Default for testing
} else {
  population <- args[1]
}

cat("=== Focal Individual Ranking Analysis ===\n")
cat(paste("Population:", population, "\n\n"))

# Define paths with updated directory structure
focal_db_dir <- file.path("output", "focal_database", population)
combined_lr_file <- file.path(focal_db_dir, paste0("focal_combined_LR_", population, ".csv"))

# We also need the unrelated database for comparison
unrelated_db_file <- file.path("output", "unrelated_database", population, paste0("unrelated_LR_", population, ".csv"))

if (!file.exists(combined_lr_file)) {
  stop(paste("Combined LR file not found:", combined_lr_file))
}

if (!file.exists(unrelated_db_file)) {
  stop(paste("Unrelated database file not found:", unrelated_db_file))
}

cat("Loading combined LR data...\n")
combined_lrs <- fread(combined_lr_file)

cat("Loading unrelated database...\n")
unrelated_lrs <- fread(unrelated_db_file)

# Check what columns we have
cat("Related pairs columns:\n")
print(names(combined_lrs))
cat("\nUnrelated database columns:\n")
print(names(unrelated_lrs))
cat("\n")

# Check what relationship types we have
cat("Available relationship types:\n")
print(unique(combined_lrs$relationship_type))
cat("\n")

# Define all populations and loci sets to analyze
populations_list <- c("AfAm", "Cauc", "Hispanic", "Asian")
loci_sets <- c("core_13", "identifiler_15", "expanded_20", "supplementary", "autosomal_29")

# Function to get the appropriate LR column for a relationship type
get_lr_column <- function(rel_type, loci_set, freq_pop) {
  if (rel_type == "full_siblings_focal") {
    return(paste0(loci_set, "_LR_full_siblings_", freq_pop))
  } else if (rel_type %in% c("half_siblings_focal", "cousins_focal", "second_cousins_focal")) {
    # Use parent_child LR for other relationships
    return(paste0(loci_set, "_LR_parent_child_", freq_pop))
  } else {
    # For any unrelated or unknown types, use parent_child as default
    return(paste0(loci_set, "_LR_parent_child_", freq_pop))
  }
}

# Function to get unrelated LR column (assuming they follow pattern: loci_set_freq_pop)
get_unrelated_lr_column <- function(loci_set, freq_pop) {
  return(paste0(loci_set, "_", freq_pop))
}

# Function to analyze individual pair rankings against unrelated database
analyze_individual_rankings <- function(related_data, unrelated_data, loci_set, freq_pop, true_pop) {
  
  cat(paste("Analyzing", loci_set, "with", freq_pop, "frequencies...\n"))
  
  # Get the unrelated LR column
  unrelated_lr_col <- get_unrelated_lr_column(loci_set, freq_pop)
  
  if (!unrelated_lr_col %in% names(unrelated_data)) {
    cat(paste("Unrelated column", unrelated_lr_col, "not found, skipping...\n"))
    return(NULL)
  }
  
  # Extract valid unrelated LR values
  unrelated_lrs_values <- unrelated_data[[unrelated_lr_col]]
  unrelated_lrs_values <- unrelated_lrs_values[is.finite(unrelated_lrs_values) & unrelated_lrs_values > 0]
  
  if (length(unrelated_lrs_values) == 0) {
    cat(paste("No valid unrelated LR values for", unrelated_lr_col, "\n"))
    return(NULL)
  }
  
  # Results storage
  individual_results <- list()
  
  # Process each related pair
  for (i in 1:nrow(related_data)) {
    rel_type <- related_data$relationship_type[i]
    pair_id <- related_data$pair_id[i]
    
    # Get the appropriate LR column for this relationship type
    related_lr_col <- get_lr_column(rel_type, loci_set, freq_pop)
    
    if (!related_lr_col %in% names(related_data)) {
      next  # Skip if column doesn't exist
    }
    
    related_lr_value <- related_data[[related_lr_col]][i]
    
    if (!is.finite(related_lr_value) || related_lr_value <= 0) {
      next  # Skip invalid values
    }
    
    # Combine this related pair LR with all unrelated LRs
    all_lrs <- c(related_lr_value, unrelated_lrs_values)
    
    # Rank all LRs (highest LR gets rank 1)
    ranks <- rank(-all_lrs, ties.method = "min")
    
    # The related pair's rank is the first one (since we put it first in the vector)
    related_pair_rank <- ranks[1]
    
    # How many unrelated pairs rank better?
    unrelated_better <- sum(unrelated_lrs_values > related_lr_value)
    
    # Store results
    individual_results[[paste(pair_id, rel_type, sep = "_")]] <- list(
      pair_id = pair_id,
      relationship_type = rel_type,
      related_lr = related_lr_value,
      rank_in_database = related_pair_rank,
      total_in_database = length(all_lrs),
      unrelated_better = unrelated_better,
      percentile = (1 - (related_pair_rank - 1) / length(all_lrs)) * 100
    )
  }
  
  if (length(individual_results) == 0) {
    cat(paste("No valid related pairs found for", loci_set, "with", freq_pop, "frequencies\n"))
    return(NULL)
  }
  
  # Convert to data frame for analysis
  results_df <- data.frame(
    pair_id = sapply(individual_results, function(x) x$pair_id),
    relationship_type = sapply(individual_results, function(x) x$relationship_type),
    related_lr = sapply(individual_results, function(x) x$related_lr),
    rank = sapply(individual_results, function(x) x$rank_in_database),
    total_database_size = sapply(individual_results, function(x) x$total_in_database),
    unrelated_better = sapply(individual_results, function(x) x$unrelated_better),
    percentile = sapply(individual_results, function(x) round(x$percentile, 2))
  )
  
  # Summary statistics
  summary_stats <- list(
    loci_set = loci_set,
    freq_population = freq_pop,
    true_population = true_pop,
    is_correct_pop = (freq_pop == true_pop),
    total_related_pairs = nrow(results_df),
    unrelated_database_size = length(unrelated_lrs_values),
    best_rank = min(results_df$rank),
    worst_rank = max(results_df$rank),
    median_rank = median(results_df$rank),
    mean_percentile = mean(results_df$percentile),
    pairs_in_top_200 = sum(results_df$rank <= 200),
    pairs_in_top_1_percent = sum(results_df$percentile >= 99),
    pairs_in_top_5_percent = sum(results_df$percentile >= 95),
    pairs_in_top_10_percent = sum(results_df$percentile >= 90),
    individual_results = results_df
  )
  
  # By relationship type summary
  by_relationship <- results_df %>%
    group_by(relationship_type) %>%
    summarise(
      count = n(),
      best_rank = min(rank),
      worst_rank = max(rank),
      median_rank = median(rank),
      mean_percentile = round(mean(percentile), 2),
      top_200 = sum(rank <= 200),
      top_200_rate = round(sum(rank <= 200) / n() * 100, 1),
      top_10_percent = sum(percentile >= 90),
      .groups = 'drop'
    )
  
  summary_stats$by_relationship = by_relationship
  
  return(summary_stats)
}

# Function to print analysis results
print_individual_analysis <- function(results) {
  if (is.null(results)) return()
  
  cat(paste("=== Individual Pair Ranking for", results$loci_set, "using", results$freq_population, "frequencies ===\n"))
  cat(paste("True population:", results$true_population, "\n"))
  cat(paste("Correct population match:", results$is_correct_pop, "\n"))
  cat(paste("Related pairs tested:", results$total_related_pairs, "\n"))
  cat(paste("Unrelated database size:", results$unrelated_database_size, "\n\n"))
  
  cat("Overall Performance:\n")
  cat(paste("Best rank achieved:", results$best_rank, "\n"))
  cat(paste("Worst rank:", results$worst_rank, "\n"))
  cat(paste("Median rank:", results$median_rank, "\n"))
  cat(paste("Mean percentile:", round(results$mean_percentile, 2), "%\n"))
  cat(paste("Pairs in top 200:", results$pairs_in_top_200, "of", results$total_related_pairs, 
            paste0("(", round(results$pairs_in_top_200/results$total_related_pairs*100, 1), "%)\n")))
  cat(paste("Pairs in top 1%:", results$pairs_in_top_1_percent, "\n"))
  cat(paste("Pairs in top 5%:", results$pairs_in_top_5_percent, "\n"))
  cat(paste("Pairs in top 10%:", results$pairs_in_top_10_percent, "\n\n"))
  
  cat("Performance by Relationship Type:\n")
  print(results$by_relationship)
  cat("\n")
  
  cat("Top 10 best-ranking pairs:\n")
  top_pairs <- results$individual_results[order(results$individual_results$rank), ][1:min(10, nrow(results$individual_results)), ]
  print(top_pairs)
  cat("\n")
}

# Run analysis for all combinations
cat("Starting individual pair ranking analysis...\n\n")

# Focus on the autosomal_29 loci set first (most comprehensive)
cat("=== AUTOSOMAL 29 LOCI ANALYSIS ===\n\n")

for (freq_pop in populations_list) {
  results <- analyze_individual_rankings(combined_lrs, unrelated_lrs, "autosomal_29", freq_pop, population)
  print_individual_analysis(results)
  cat(paste(rep("=", 80), collapse = ""), "\n\n")
}

# Then check other loci sets using the correct population frequencies
cat("=== LOCI SET COMPARISON (using correct population frequencies) ===\n\n")

for (loci_set in loci_sets) {
  if (loci_set == "autosomal_29") next  # Already done above
  
  results <- analyze_individual_rankings(combined_lrs, unrelated_lrs, loci_set, population, population)
  print_individual_analysis(results)
  cat(paste(rep("=", 80), collapse = ""), "\n\n")
}

# Summary comparison across all loci sets and frequency populations
cat("=== SUMMARY COMPARISON ===\n\n")

summary_results <- list()

for (loci_set in loci_sets) {
  for (freq_pop in populations_list) {
    results <- analyze_individual_rankings(combined_lrs, unrelated_lrs, loci_set, freq_pop, population)
    
    if (!is.null(results)) {
      summary_key <- paste(loci_set, freq_pop, sep = "_")
      summary_results[[summary_key]] <- list(
        loci_set = loci_set,
        freq_population = freq_pop,
        is_correct_pop = results$is_correct_pop,
        best_rank = results$best_rank,
        median_rank = results$median_rank,
        mean_percentile = results$mean_percentile,
        total_pairs = results$total_related_pairs,
        top_200 = results$pairs_in_top_200,
        top_200_rate = round(results$pairs_in_top_200 / results$total_related_pairs * 100, 1),
        top_10_percent = results$pairs_in_top_10_percent,
        top_5_percent = results$pairs_in_top_5_percent,
        top_1_percent = results$pairs_in_top_1_percent
      )
    }
  }
}

if (length(summary_results) > 0) {
  summary_df <- data.frame(
    loci_set = sapply(summary_results, function(x) x$loci_set),
    freq_pop = sapply(summary_results, function(x) x$freq_population),
    correct_pop = sapply(summary_results, function(x) x$is_correct_pop),
    best_rank = sapply(summary_results, function(x) x$best_rank),
    median_rank = sapply(summary_results, function(x) x$median_rank),
    mean_percentile = sapply(summary_results, function(x) round(x$mean_percentile, 1)),
    total_pairs = sapply(summary_results, function(x) x$total_pairs),
    top_200 = sapply(summary_results, function(x) x$top_200),
    top_200_rate = sapply(summary_results, function(x) paste0(x$top_200_rate, "%")),
    top_10_pct = sapply(summary_results, function(x) x$top_10_percent),
    top_5_pct = sapply(summary_results, function(x) x$top_5_percent),
    top_1_pct = sapply(summary_results, function(x) x$top_1_percent)
  )
  
  # Sort by best rank, then by mean percentile
  summary_df <- summary_df[order(summary_df$best_rank, -summary_df$mean_percentile), ]
  
  cat("Summary of individual pair ranking performance:\n")
  print(summary_df)
  
  # Save summary to file
  output_file <- file.path(focal_db_dir, paste0("focal_individual_ranking_summary_", population, ".csv"))
  write.csv(summary_df, output_file, row.names = FALSE)
  cat(paste("\nDetailed summary saved to:", output_file, "\n"))
  
  # Also save detailed results for the best performing combination
  best_combo <- summary_df[1, ]
  cat(paste("\nBest performing combination:", best_combo$loci_set, "with", best_combo$freq_pop, "frequencies\n"))
  
} else {
  cat("No valid results found for summary\n")
}

cat("\n=== Analysis Complete ===\n")