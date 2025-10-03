#!/usr/bin/env Rscript
# analyze_lr_outputs.R
# Process chunked LR and combined_LR outputs for plotting
# Fixed version that handles multiple tested_relationships per pair

library(tidyverse)
library(data.table)

# Set up output directory
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  output_dir <- args[1]
} else {
  output_dir <- paste0("output/lr_analysis_", format(Sys.time(), "%Y%m%d_%H%M%S"))
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Output directory:", output_dir, "\n")

# Define populations and relationships
populations <- c("AfAm", "Cauc", "Hispanic", "Asian")
relationships <- c("parent_child", "full_siblings", "half_siblings", 
                  "cousins", "second_cousins", "unrelated")

# Read all LR files
cat("Reading LR files...\n")
all_lr <- list()
for (pop in populations) {
  for (rel in relationships) {
    files <- list.files("output/LR/", 
                       pattern = paste0("^LR_", pop, "_", rel, "_n1000_chunk.*\\.csv$"),
                       full.names = TRUE)
    if (length(files) > 0) {
      cat("Reading", length(files), "files for", pop, rel, "\n")
      temp_list <- list()
      for (f in files) {
        temp_list[[f]] <- fread(f)
      }
      all_lr[[paste(pop, rel, sep="_")]] <- rbindlist(temp_list)
      rm(temp_list)
    }
  }
}
all_lr <- rbindlist(all_lr)
gc()

# Read all combined LR files
cat("Reading combined LR files...\n")
all_combined <- list()
for (pop in populations) {
  for (rel in relationships) {
    files <- list.files("output/combined_LR/", 
                       pattern = paste0("^combined_LR_", pop, "_", rel, "_n1000_chunk.*\\.csv$"),
                       full.names = TRUE)
    if (length(files) > 0) {
      cat("Reading", length(files), "files for", pop, rel, "\n")
      temp_list <- list()
      for (f in files) {
        dt <- fread(f)
        # Ensure combined_LR is numeric immediately after reading
        if ("combined_LR" %in% names(dt)) {
          dt[, combined_LR := as.numeric(combined_LR)]
        }
        temp_list[[f]] <- dt
      }
      all_combined[[paste(pop, rel, sep="_")]] <- rbindlist(temp_list)
      rm(temp_list)
    }
  }
}
all_combined <- rbindlist(all_combined)
gc()

# Add is_correct_pop column (keeping original column names)
all_lr[, is_correct_pop := (population == tested_population)]
all_combined[, is_correct_pop := (population == tested_population)]

# Process each population
for (pop in populations) {
  cat("\nProcessing", pop, "...\n")
  
  # Create output directory for this population
  pop_dir <- file.path(output_dir, paste0(pop, "_summary"))
  dir.create(pop_dir, showWarnings = FALSE)
  
  # Filter data for this population
  pop_combined <- all_combined[population == pop]
  
  if (nrow(pop_combined) == 0) {
    cat("No data for", pop, "\n")
    next
  }
  
  # Save combined LR data (this is what plots need)
  fwrite(pop_combined, file.path(pop_dir, paste0("sim_summary_genotypes_", pop, ".csv")))
  
  # Calculate summary statistics - NOW INCLUDING tested_relationship as grouping variable
  summary_stats <- pop_combined[, .(
    n = .N,
    mean_LR = mean(combined_LR, na.rm = TRUE),
    median_LR = median(combined_LR, na.rm = TRUE),
    sd_LR = sd(combined_LR, na.rm = TRUE),
    min_LR = min(combined_LR, na.rm = TRUE),
    max_LR = max(combined_LR, na.rm = TRUE),
    lower_95 = quantile(combined_LR, 0.025, na.rm = TRUE),
    upper_95 = quantile(combined_LR, 0.975, na.rm = TRUE)
  ), by = .(known_relationship, population, loci_set, tested_relationship, tested_population, is_correct_pop)]
  
  fwrite(summary_stats, file.path(pop_dir, paste0("sim_lr_summary_stats_", pop, ".csv")))
  
  # Calculate ratio statistics (wrong/correct) - NOW INCLUDING tested_relationship
  # Create unique pair identifier to avoid cartesian join
  pop_combined[, unique_pair_id := paste(batch_id, pair_id, sep="_")]
  
  correct_lr <- pop_combined[is_correct_pop == TRUE, .(
    population, known_relationship, unique_pair_id, loci_set, tested_relationship, correct_LR = combined_LR
  )]
  
  wrong_lr <- pop_combined[is_correct_pop == FALSE, .(
    population, known_relationship, unique_pair_id, loci_set, tested_relationship, tested_population, wrong_LR = combined_LR
  )]
  
  # Now merge will work because tested_relationship makes keys unique
  ratio_data <- merge(wrong_lr, correct_lr, 
                     by = c("population", "known_relationship", "unique_pair_id", "loci_set", "tested_relationship"))
  ratio_data[, ratio := wrong_LR / correct_LR]
  
  ratio_summary <- ratio_data[, .(
    n = .N,
    mean_ratio = mean(ratio, na.rm = TRUE),
    median_ratio = median(ratio, na.rm = TRUE),
    sd_ratio = sd(ratio, na.rm = TRUE),
    min_ratio = min(ratio, na.rm = TRUE),
    max_ratio = max(ratio, na.rm = TRUE),
    lower_95 = quantile(ratio, 0.025, na.rm = TRUE),
    upper_95 = quantile(ratio, 0.975, na.rm = TRUE)
  ), by = .(known_relationship, population, loci_set, tested_relationship, tested_population)]
  
  fwrite(ratio_summary, file.path(pop_dir, paste0("sim_lr_ratio_summary_", pop, ".csv")))
  
  # Calculate cutoffs based on unrelated pairs - NOW INCLUDING tested_relationship
  unrelated <- pop_combined[known_relationship == "unrelated"]
  
  cutoffs <- unrelated[, .(
    fixed_cutoff = 1.0,
    cutoff_1 = quantile(combined_LR, 0.99, na.rm = TRUE),
    cutoff_0_1 = quantile(combined_LR, 0.999, na.rm = TRUE),
    cutoff_0_01 = quantile(combined_LR, 0.9999, na.rm = TRUE),
    n_unrelated = .N
  ), by = .(loci_set, tested_relationship, tested_population)]
  
  fwrite(cutoffs, file.path(pop_dir, paste0("sim_cutoffs_", pop, ".csv")))
  
  # Calculate proportions exceeding cutoffs - NOW INCLUDING tested_relationship
  with_cutoffs <- merge(pop_combined, cutoffs, 
                        by = c("loci_set", "tested_relationship", "tested_population"),
                        all.x = TRUE)
  
  proportions <- with_cutoffs[known_relationship != "unrelated", .(
    n_pairs = .N,
    proportion_exceeding_fixed = mean(combined_LR > fixed_cutoff, na.rm = TRUE),
    proportion_exceeding_1 = mean(combined_LR > cutoff_1, na.rm = TRUE),
    proportion_exceeding_0_1 = mean(combined_LR > cutoff_0_1, na.rm = TRUE),
    proportion_exceeding_0_01 = mean(combined_LR > cutoff_0_01, na.rm = TRUE)
  ), by = .(population, known_relationship, loci_set, tested_relationship, tested_population, is_correct_pop)]
  
  fwrite(proportions, file.path(pop_dir, paste0("sim_proportions_exceeding_cutoffs_", pop, ".csv")))
  
  cat("Saved files for", pop, "\n")
}

# Quick summary
cat("\n=== Summary ===\n")
cat("Total single-locus LRs:", nrow(all_lr), "\n")
cat("Total combined LRs:", nrow(all_combined), "\n")
cat("Unique pairs:", length(unique(all_combined$pair_id)), "\n")

# Sample size table - NOW INCLUDING tested_relationship
sample_sizes <- all_combined[, .(n = uniqueN(pair_id)), by = .(population, known_relationship, tested_relationship)]
sample_table <- dcast(sample_sizes, known_relationship + tested_relationship ~ population, value.var = "n")
cat("\nSample sizes by population and tested relationship:\n")
print(sample_table)

# Show relationship combinations
rel_combinations <- all_combined[, .(
  n_pairs = uniqueN(pair_id),
  n_rows = .N
), by = .(known_relationship, tested_relationship)]
cat("\nKnown vs tested relationship combinations:\n")
print(rel_combinations)

cat("\nAnalysis complete. Output saved to:", output_dir, "\n")
cat("\nTo generate plots, run:\n")
cat("# For matched relationship analysis (tested == known):\n")
cat("Rscript code/plots_known_NEA.R", output_dir, "\n")
cat("# For mismatched population analysis:\n")
cat("Rscript code/plots_known_vs_tested_NEA.R", output_dir, "\n")
cat("\nNote: You may need to modify your plotting scripts to handle the new tested_relationship column\n")
cat("or filter the data to specific scenarios (e.g., tested_relationship == known_relationship)\n")