#!/usr/bin/env Rscript

# Load Required Libraries
suppressMessages(suppressWarnings({
  library(tidyverse)  # Includes ggplot2, dplyr, tidyr, etc.
  library(data.table) # Required for efficient fread()
  library(scales)     # For number formatting
}))

# Helper function for logging
log_message <- function(message) {
  cat(paste0("[", Sys.time(), "] ", message, "\n"))
}

# --- Argument Parsing and Setup ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript code/plots_matched.R <input_dir> [output_dir]")
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

# Define relationship types, colors, and labels
relationship_order <- c("parent_child", "full_siblings", "half_siblings", 
                        "cousins", "second_cousins", "unrelated")
relationship_labels <- c("Parent-Child", "Full Siblings", "Half Siblings", 
                         "Cousins", "Second Cousins", "Unrelated")
names(relationship_labels) <- relationship_order

loci_set_order <- c("core_13", "identifiler_15", "expanded_20", "supplementary", "autosomal_29")

# Define color palette for Loci Sets
loci_colors <- c(
  "core_13" = "#1b9e77",        # Teal/Green
  "identifiler_15" = "#d95f02", # Orange
  "expanded_20" = "#7570b3",    # Purple
  "supplementary" = "#e7298a",  # Pink
  "autosomal_29" = "#66a61e"    # Light Green
)

# --- Data Loading Function ---

load_combined_data <- function(input_dir) {
  
  # --- 1. Load Raw Combined LR Data (Strictly Matched) ---
  raw_lrs_path <- file.path(input_dir, COMBINED_LR_MATCH_FILE)
  combined_lrs_match <- tryCatch({
    dt <- fread(raw_lrs_path)
    
    # The raw file only contains matched data (pop and rel), but re-filter for safety
    dt %>% 
      as_tibble() %>%
      filter(is_correct_pop == TRUE) %>% filter(known_relationship == tested_relationship) %>%
      mutate(
        known_relationship = factor(known_relationship, levels = relationship_order),
        loci_set = factor(loci_set, levels = loci_set_order),
        combined_LR = as.numeric(combined_LR)
      )
    
  }, error = function(e) {
    log_message(paste("Could not load raw matched LR data:", e$message))
    return(NULL)
  })
  
  # --- 2. Load Summary Stats ---
  summary_stats_path <- file.path(input_dir, SUMMARY_STATS_FILE)
  summary_stats <- tryCatch({
    dt <- fread(summary_stats_path)
    
    # Convert character-based numeric columns now
    dt <- dt %>%
      mutate(across(c(mean_LR, median_LR, sd_LR, min_LR, max_LR, lower_95, upper_95), as.numeric))
    
    # CRITICAL FILTER: Keep only Population-Matched AND Relationship-Matched (known == tested)
    dt_matched <- dt %>% 
      filter(is_correct_pop == TRUE) %>% 
      filter(known_relationship == tested_relationship)
    
    dt_matched %>% 
      as_tibble() %>%
      mutate(
        known_relationship = factor(known_relationship, levels = relationship_order),
        loci_set = factor(loci_set, levels = loci_set_order),
      )
    
  }, error = function(e) {
    log_message(paste("Could not load summary stats:", e$message))
    return(NULL)
  })
  
  # --- 3. Load Proportions Exceeding Cutoffs ---
  proportions_path <- file.path(input_dir, PROPORTIONS_FILE)
  proportions <- tryCatch({
    dt <- fread(proportions_path)
    
    # CRITICAL FILTER: Keep only Population-Matched AND Relationship-Matched (known == tested)
    dt_matched <- dt %>% 
      filter(is_correct_pop == TRUE) %>% 
      filter(known_relationship == tested_relationship)
    
    dt_matched %>% 
      as_tibble() %>%
      mutate(
        known_relationship = factor(known_relationship, levels = relationship_order),
        loci_set = factor(loci_set, levels = loci_set_order),
        population = factor(population, levels = c("AfAm", "Asian", "Cauc", "Hispanic"))
      )
    
  }, error = function(e) {
    log_message(paste("Could not load proportions data:", e$message))
    return(NULL)
  })
  
  # --- 4. Load Cutoffs ---
  cutoffs_path <- file.path(input_dir, CUTOFFS_FILE)
  cutoffs <- tryCatch({
    dt <- fread(cutoffs_path)
    dt %>% as_tibble()
  }, error = function(e) {
    log_message(paste("Could not load cutoffs data:", e$message))
    return(NULL)
  })
  
  population_relationship_tallies <- combined_lrs_match %>% group_by(population, known_relationship) %>% tally()
  
  return(list(
    combined_lrs_match = combined_lrs_match,
    summary_stats = summary_stats,
    proportions = proportions,
    cutoffs = cutoffs,
    population_relationship_tallies = population_relationship_tallies
  ))
}

# --- Plotting Functions (Transferred and Adapted) ---

plot_lr_distributions <- function(lrs_data, population_relationship_tallies) {
  # Calculate total samples per population from the tallies
  population_totals <- population_relationship_tallies %>%
    group_by(population) %>%
    summarize(total_n = sum(n), .groups = 'drop')
  
  # Create legend labels with population names and total sample counts
  legend_labels <- setNames(
    paste0(population_totals$population, " (n=", format(population_totals$total_n, big.mark = ","), ")"),
    population_totals$population
  )
  
  # Define a color palette
  pop_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
  names(pop_colors) <- levels(lrs_data$population)
  
  plot_obj <- ggplot(lrs_data, aes(x = known_relationship, y = combined_LR, fill = population)) +
    geom_boxplot(position = position_dodge(width = 0.85), alpha = 0.7, outlier.size = 0.5) +
    facet_wrap(~ loci_set, scales = "fixed") +
    scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_fill_manual(values = pop_colors, labels = legend_labels) +
    labs(
      title = "LR Distributions Across Populations and Relationship Types",
      subtitle = "Using Population-Matched Allele Frequencies",
      x = "Relationship Type",
      y = "LR (log scale)",
      fill = "Population"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.title = element_text(size = 14),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.position = "bottom",
      strip.text = element_text(size = 12)
    )
  return(plot_obj)
}

plot_mean_lr <- function(summary_data, population_relationship_tallies) {
  # Calculate total samples per population from the tallies
  population_totals <- population_relationship_tallies %>%
    group_by(population) %>%
    summarize(total_n = sum(n), .groups = 'drop')
  
  # Create legend labels with population names and total sample counts
  legend_labels <- setNames(
    paste0(population_totals$population, " (n=", format(population_totals$total_n, big.mark = ","), ")"),
    population_totals$population
  )
  
  # Calculate total samples per relationship type for facet labels
  relationship_totals <- population_relationship_tallies %>%
    group_by(known_relationship) %>%
    summarize(total_n = sum(n), .groups = 'drop')
  
  # Create facet labels with relationship types and total sample counts
  facet_labels <- setNames(
    paste0(relationship_totals$known_relationship, "\n(n=", format(relationship_totals$total_n, big.mark = ","), ")"),
    relationship_totals$known_relationship
  )
  
  # Define a color palette
  pop_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
  names(pop_colors) <- levels(summary_data$population)
  
  plot_obj <- ggplot(summary_data, aes(x = loci_set, y = mean_LR, group = population, color = population)) +
      geom_line(size = 1.2) +
      geom_point(size = 3) +
      facet_wrap(~ known_relationship, scales = "fixed", ncol = 2,
                 labeller = labeller(known_relationship = facet_labels)) +
      scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      scale_color_manual(values = pop_colors, labels = legend_labels) +
      labs(
        title = "Mean LR Across Populations and Relationship Types",
        subtitle = "Using Population-Matched Allele Frequencies",
        x = "Loci Set",
        y = "Mean Combined LR (log scale)",
        color = "Population"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.position = "bottom",
        strip.text = element_text(size = 12)
      )
    
  return(plot_obj)
}

plot_proportions_exceeding_cutoffs <- function(proportions_data) {
  # Define a color palette
  pop_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
  names(pop_colors) <- levels(proportions_data$population)
  
  # Convert to long format for easy facetting
    proportions_long <- proportions_data %>%
      pivot_longer(cols = starts_with("proportion_exceeding"),
                   names_to = "Cutoff_Type", values_to = "Proportion",
                   names_prefix = "proportion_exceeding_")
    
    proportions_long$Cutoff_Type <- factor(
      proportions_long$Cutoff_Type, 
      levels = c("fixed", "1", "0_1", "0_01"),
      labels = c("Fixed Cutoff (1.00)", "1% FPR", "0.1% FPR", "0.01% FPR")
    )
    
    plot_obj <- ggplot(proportions_long, aes(x = known_relationship, y = Proportion, fill = population)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
      facet_grid(loci_set ~ Cutoff_Type) +
      scale_fill_manual(values = pop_colors) +
      labs(
        title = "Proportions Exceeding Likelihood Cut-offs",
        subtitle = "Using Population-Matched Allele Frequencies",
        x = "Relationship Type",
        y = "Proportion Exceeding Cut-off",
        fill = "Population"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.position = "bottom",
        strip.text = element_text(size = 9)
      )
  return(plot_obj)
}

plot_heatmap_proportions <- function(proportions_data) {
  
  # Convert the proportions to long format
  proportions_long <- proportions_data %>%
    pivot_longer(cols = starts_with("proportion_exceeding"),
                 names_to = "Cutoff_Type", values_to = "Proportion",
                 names_prefix = "proportion_exceeding_")
  
  proportions_long$Cutoff_Type <- factor(
    proportions_long$Cutoff_Type, 
    levels = c("fixed", "1", "0_1", "0_01"),
    labels = c("Fixed Cutoff (1.00)", "1% FPR", "0.1% FPR", "0.01% FPR")
  )
  
  # Create a custom color gradient from blue to red
  color_gradient <- colorRampPalette(c("blue", "white", "red"))(100)
  
  plot_obj <- ggplot(proportions_long, aes(x = loci_set, y = known_relationship, fill = Proportion)) +
    geom_tile() +
    facet_grid(population ~ Cutoff_Type) +
    scale_fill_gradientn(colors = color_gradient, limits = c(0, 1)) +
    labs(
      title = "Heat Map of Proportions Exceeding Likelihood Cut-offs",
      subtitle = "Using Population-Matched Allele Frequencies",
      x = "Loci Set",
      y = "Relationship Type",
      fill = "Proportion"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.title = element_text(size = 14),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.position = "bottom",
      strip.text = element_text(size = 10)
    )
  
  return(plot_obj)
}


# --- Plot Saving Utility ---
save_plot <- function(plot_obj, filename, width = 10, height = 8) {
  if (!is.null(plot_obj)) {
    ggsave(file.path(output_dir, filename), plot = plot_obj, width = width, height = height, bg = "white")
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
  
  # --- 1. Raw LR Distributions (Box/Violin) ---
  if (!is.null(all_data$combined_lrs_match) && nrow(all_data$combined_lrs_match) > 0) {
    log_message("Creating LR distribution plots...")
    
    # Box plot
    lr_boxplot <- plot_lr_distributions(all_data$combined_lrs_match, all_data$population_relationship_tallies)
    save_plot(lr_boxplot, "lr_distributions_boxplot_matched.png", width = 12, height = 8)
    
  } else {
    log_message("No combined_lrs_match data available for distribution plots.")
  }
  
  # --- 2. Mean LR Plot ---
  if (!is.null(all_data$summary_stats) && nrow(all_data$summary_stats) > 0) {
    log_message("Creating Mean LR plot...")
    
    mean_lr_plot <- plot_mean_lr(all_data$summary_stats, all_data$population_relationship_tallies)
    save_plot(mean_lr_plot, "mean_combined_lr_matched.png", width = 12, height = 8)
  } else {
    log_message("No summary_stats data available for Mean LR plot.")
  }
  
  # --- 3. Proportions Plots ---
  if (!is.null(all_data$proportions) && nrow(all_data$proportions) > 0) {
    log_message("Creating Proportions Exceeding Cutoffs plots...")
    
    # Proportions bar plot
    proportions_plot <- plot_proportions_exceeding_cutoffs(all_data$proportions)
    save_plot(proportions_plot, "proportions_exceeding_cutoffs_matched.png", width = 10, height = 12)
    
    # Proportions heatmap
    heatmap_plot <- plot_heatmap_proportions(all_data$proportions)
    save_plot(heatmap_plot, "heatmap_proportions_fixed_cutoff_matched.png", width = 8, height = 8)
  } else {
    log_message("No proportions data available for plotting.")
  }
  
  log_message("All available plots created and saved successfully.")
}

log_message("Script finished.")