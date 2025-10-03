#!/usr/bin/env Rscript

# Load Required Libraries
suppressMessages(suppressWarnings({
  library(tidyverse)  # Includes ggplot2, dplyr, tidyr, etc.
  library(data.table) # Required for efficient fread()
}))

# Helper function for logging
log_message <- function(message) {
  cat(paste0("[", Sys.time(), "] ", message, "\n"))
}

# --- Argument Parsing and Setup ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript code/plots_known_NEA.R <input_dir> [output_dir]")
}

input_dir <- args[1]
log_message(paste("Input directory:", input_dir))

if (length(args) >= 2) {
  output_subdir <- args[2]
  output_dir <- file.path("output", output_subdir)
} else {
  timestamp <- format(Sys.time(), "%Y%m%d")
  output_dir <- file.path(input_dir, paste0("matched_pop_plots_", timestamp)) 
}
log_message(paste("Output directory:", output_dir))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Define file names
COMBINED_LR_MATCH_FILE <- "combined_LR_match.csv"
SUMMARY_STATS_FILE <- "combined_LR_summary_stats.csv"
PROPORTIONS_FILE <- "combined_LR_exceeding_cutoffs.csv"
CUTOFFS_FILE <- "combined_LR_cutoffs.csv" 

# Define relationship types and ensure correct ordering for plotting
relationship_order <- c("parent_child", "full_siblings", "half_siblings", "cousins", "second_cousins", "unrelated")
relationship_labels <- c("Parent-Child", "Full Siblings", "Half Siblings", "Cousins", "Second Cousins", "Unrelated")

# Define loci sets and ensure correct ordering for plotting
loci_set_order <- c("core_13", "identifiler_15", "expanded_20", "supplementary", "autosomal_29")
loci_set_labels <- c("Core 13", "Identifiler 15", "Expanded 20", "Supplementary", "Autosomal 29")
# ---------------------------------------------------------------------


# --- Data Loading Function ---
load_combined_data <- function(input_dir) {
  
  # --- 1. Load Raw Combined LR Data (Matched Only) ---
  raw_lrs_path <- file.path(input_dir, COMBINED_LR_MATCH_FILE)
  combined_lrs_match <- tryCatch({
    dt <- fread(raw_lrs_path)
    
    # Ensure correct column types and order (This file only contains matched data)
    dt %>% 
      as_tibble() %>%
      filter(is_correct_pop == TRUE) %>% filter(known_relationship == tested_relationship) %>% 
      mutate(
        known_relationship = factor(known_relationship, levels = relationship_order),
        loci_set = factor(loci_set, levels = loci_set_order)
      )
      
  }, error = function(e) {
    log_message(paste("Could not load raw matched LR data:", e$message))
    return(NULL)
  })
  
  # --- 2. Load Summary Stats ---
  summary_stats_path <- file.path(input_dir, SUMMARY_STATS_FILE)
  summary_stats <- tryCatch({
    dt <- fread(summary_stats_path)
    
    # Note: Summary stats file is filtered for is_correct_pop == TRUE during creation, 
    # but we add a safety filter here in case the input file isn't perfectly clean.
    dt_matched <- dt %>% filter(is_correct_pop == TRUE) %>% filter(known_relationship == tested_relationship)
    
    dt_matched %>% 
      as_tibble() %>%
      mutate(
        known_relationship = factor(known_relationship, levels = relationship_order),
        loci_set = factor(loci_set, levels = loci_set_order)
      )
      
  }, error = function(e) {
    log_message(paste("Could not load summary stats:", e$message))
    return(NULL)
  })
  
  # --- 3. Load Proportions Exceeding Cutoffs ---
  proportions_path <- file.path(input_dir, PROPORTIONS_FILE)
  proportions <- tryCatch({
    dt <- fread(proportions_path)
    
    # Safety filter for matched population
    dt_matched <- dt %>% filter(is_correct_pop == TRUE) %>% filter(known_relationship == tested_relationship)
    
    dt_matched %>% 
      as_tibble() %>%
      mutate(
        known_relationship = factor(known_relationship, levels = relationship_order),
        loci_set = factor(loci_set, levels = loci_set_order)
      )
      
  }, error = function(e) {
    log_message(paste("Could not load proportions:", e$message))
    return(NULL)
  })
  
  # --- 4. Load Cutoffs ---
  cutoffs_path <- file.path(input_dir, CUTOFFS_FILE)
  cutoffs <- tryCatch({
    dt <- fread(cutoffs_path)
    
    dt %>% 
      as_tibble() %>%
      mutate(
        tested_relationship = factor(tested_relationship, levels = relationship_order),
        loci_set = factor(loci_set, levels = loci_set_order)
      )
      
  }, error = function(e) {
    log_message(paste("Could not load cutoffs:", e$message))
    return(NULL)
  })
  
  return(list(
    combined_lrs_match = combined_lrs_match, # <-- ADDED RAW LR DATA
    summary_stats = summary_stats,
    proportions = proportions,
    cutoffs = cutoffs 
  ))
}
# ---------------------------------------------------


# --- Plotting Functions Placeholder ---
save_plot <- function(plot_obj, filename, width=10, height=8) {
  if (!is.null(plot_obj)) {
    ggsave(file.path(output_dir, filename), plot = plot_obj, width = width, height = height)
    log_message(paste("Saved plot:", filename))
  }
}
# -------------------------------------------------------------------------------


# Main execution
log_message("Starting matched population plotting process...")

all_data <- tryCatch({
  load_combined_data(input_dir)
}, error = function(e) {
  log_message(paste("Critical error during data loading:", e$message))
  return(NULL)
})


if (!is.null(all_data)) {
  
  if (!is.null(all_data$combined_lrs_match) && nrow(all_data$combined_lrs_match) > 0) {
    log_message("Raw matched LR data loaded successfully for plotting LR distributions (e.g., box plots).")
    # Example plotting function calls using the raw data:
    # lr_distribution_plot <- plot_lr_distributions(all_data$combined_lrs_match, all_data$cutoffs)
    # save_plot(lr_distribution_plot, "lr_distributions_matched.png")
  }
  
  if (!is.null(all_data$summary_stats) && nrow(all_data$summary_stats) > 0) {
    log_message("Summary statistics data loaded successfully.")
  }
  
  if (!is.null(all_data$proportions) && nrow(all_data$proportions) > 0) {
    log_message("Proportions data loaded successfully.")
  }
  
  if (!is.null(all_data$cutoffs) && nrow(all_data$cutoffs) > 0) {
    log_message("Cutoffs data loaded successfully.")
  }
  
  log_message("Script execution finished. Review output directory for saved plots.")
} else {
  log_message("Script terminated due to critical data loading error.")
}