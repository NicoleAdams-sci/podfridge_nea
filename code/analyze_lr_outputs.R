#!/usr/bin/env Rscript
# analyze_lr_outputs.R
# Process chunked LR and combined_LR outputs for plotting
# Fixed version that handles multiple tested_relationships per pair

library(tidyverse)
library(data.table)

# Source required modules
source("code/module9_combinedLR_stats_functions.R")

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
populations <- c("AfAm", "Cauc", "Hispanic", "Asian", "all")
relationships <- c("parent_child", "full_siblings", "half_siblings", 
                  "cousins", "second_cousins", "unrelated")

# Read in combined_LR files
all_combined <- list()
for (pop in populations) {
  for (rel in relationships) {
    files <- list.files("output/combined_LR", 
                        pattern = paste0("^combined_LR_", pop, "_", rel, "_n1000_chunk.*\\.csv$"),
                        full.names = TRUE)
    if (length(files) > 0) {
      cat("Reading", length(files), "files for", pop, rel, "\n")
      temp_list <- list()
      for (f in files) {
        dt <- fread(f)
        temp_list[[f]] <- dt
      }
      all_combined[[paste(pop, rel, sep="_")]] <- rbindlist(temp_list)
      rm(temp_list)
    }
  }
}
all_combined <- rbindlist(all_combined)
gc()


# make sure combined_LR is numeric
all_combined <- all_combined %>% mutate(combined_LR = as.numeric(combined_LR))

# save all_combined
all_combined_file <- file.path(output_dir, "combined_LR_all.rds")
saveRDS(all_combined, all_combined_file, compress = "gzip")

# 1. Define the data frame for Population Match ONLY (used for summary stats and proportions)
# This data frame includes ALL tested relationships, as long as the population is correct.
all_combined_pop_match <- all_combined %>% filter(is_correct_pop == TRUE)

# 2. Define the data frame for STRICT Match (used for the raw LR distribution plots)
all_combined_strict_match <- all_combined_pop_match %>% filter(known_relationship == tested_relationship) %>% mutate(unique_pair_id = str_c(batch_id, pair_id, sep = "_"))

# 3. Define the data frame for Population MISMATCH (used for ratio analysis), relationship is known
all_combined_mismatch <- all_combined %>% filter(is_correct_pop == FALSE)

# 4. Save the combined LR files (Uses the STRICT match data)
all_combined_match_file <- file.path(output_dir, "combined_LR_match.csv")
fwrite(all_combined_strict_match, all_combined_match_file) # <-- Use the STRICT match data here

# Clean up memory
rm(all_combined_strict_match)
gc()

all_combined_mismatch_file <- file.path(output_dir, "combined_LR_mismatch.csv")
fwrite(all_combined_mismatch, all_combined_mismatch_file)

# Clean up memory
rm(all_combined_mismatch)
gc()

# apply functions from module 9 - analyses on combined_LR data
# Use the POPULATION-MATCH ONLY data for summary stats
summary_stats <- calculate_summary_stats(all_combined_pop_match)
summary_stats_file <- file.path(output_dir, "combined_LR_summary_stats.csv")
fwrite(summary_stats, summary_stats_file)

# Clean up memory
rm(summary_stats)
gc()

# Uses FULL combined data
ratios <- calculate_ratio_stats(all_combined)
ratio_summary_file <- file.path(output_dir, "combined_LR_ratio_summary.csv")
fwrite(ratios$ratio_summary, ratio_summary_file)

ratios_raw_file <- file.path(output_dir, "combined_LR_ratios_raw.csv")
fwrite(ratios$combined_lrs_ratio, ratios_raw_file)

# Clean up memory
rm(ratios)
gc()