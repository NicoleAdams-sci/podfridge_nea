#!/usr/bin/env Rscript

# Load Required Libraries
suppressMessages(suppressWarnings({
  library(tidyverse)  # Includes ggplot2, dplyr, tidyr, etc.
}))

# Helper function for logging
log_message <- function(message) {
  cat(paste0("[", Sys.time(), "] ", message, "\n"))
}

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript code/plots_known_NEA.R <input_dir> [output_dir]")
}

input_dir <- args[1]
log_message(paste("Input directory:", input_dir))

# If output directory is provided, use it, otherwise create one based on the input directory
if (length(args) >= 2) {
  output_subdir <- args[2]
  output_dir <- file.path("output", output_subdir)
} else {
  timestamp <- format(Sys.time(), "%Y%m%d")
  output_dir <- file.path("output", paste0("matched_pop_plots_", timestamp))
}
log_message(paste("Output directory:", output_dir))

# Create output directory if it doesn't exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Define populations
populations <- c("AfAm", "Cauc", "Hispanic", "Asian")

# Define relationship types and ensure correct ordering for plotting
relationship_order <- c("parent_child", "full_siblings", "half_siblings", "cousins", "second_cousins", "unrelated")
relationship_labels <- c("Parent-Child", "Full Siblings", "Half Siblings", "Cousins", "Second Cousins", "Unrelated")

# Define loci sets and ensure correct ordering for plotting
loci_set_order <- c("core_13", "identifiler_15", "expanded_20", "supplementary", "autosomal_29")
loci_set_labels <- c("Core 13", "Identifiler 15", "Expanded 20", "Supplementary", "Autosomal 29")


# Load data from each population - BUT FILTER ONLY FOR MATCHED POPULATIONS
load_combined_data <- function() {
  log_message("Loading combined LR data for all populations...")
  
  # Initialize empty data frames
  all_combined_lrs <- data.frame()
  all_cutoffs <- data.frame()
  all_proportions <- data.frame()
  all_summary_stats <- data.frame()
  
  # Track how many files were successfully loaded
  files_loaded <- 0
  
  for (pop in populations) {
    log_message(paste("Processing population:", pop))
    
    # Define the population-specific directory
    pop_dir <- file.path(input_dir, paste0(pop, "_summary"))
    
    
    # Load combined LRs
    combined_lrs_file <- file.path(pop_dir, paste0("sim_summary_genotypes_", pop, ".csv"))
    data <- read.csv(combined_lrs_file)
    
    if (!is.null(data)) {
      # ONLY keep rows where population matches frequency population (is_correct_pop == TRUE)
      data <- data %>% filter(is_correct_pop == TRUE)
      
      if (nrow(data) > 0) {
        files_loaded <- files_loaded + 1
        all_combined_lrs <- bind_rows(all_combined_lrs, data)
      }
    }
    
    # Load cutoffs for the correct population
    cutoffs_file <- file.path(pop_dir, paste0("sim_cutoffs_", pop, ".csv"))
    data <- read.csv(cutoffs_file)
    
    if (!is.null(data)) {
      # Filter cutoffs to only include the matching population
      data <- data %>% filter(freq_population == pop)
      all_cutoffs <- bind_rows(all_cutoffs, data)
    }
    
    # Load proportions
    proportions_file <- file.path(pop_dir, paste0("sim_proportions_exceeding_cutoffs_", pop, ".csv"))
    data <- read.csv(proportions_file)
    
    if (!is.null(data)) {
      # ONLY keep rows where population matches frequency population (is_correct_pop == TRUE)
      data <- data %>% filter(is_correct_pop == TRUE)
      all_proportions <- bind_rows(all_proportions, data)
    }
    
    # Load summary stats
    summary_stats_file <- file.path(pop_dir, paste0("sim_lr_summary_stats_", pop, ".csv"))
    data <- read.csv(summary_stats_file)
    
    if (!is.null(data)) {
      # ONLY keep rows where population matches frequency population (is_correct_pop == TRUE)
      data <- data %>% filter(is_correct_pop == TRUE)
      all_summary_stats <- bind_rows(all_summary_stats, data)
    }
  }
  
  # Check if we have data
  if (nrow(all_combined_lrs) == 0) {
    log_message("ERROR: No combined LR data was loaded. Please check the input directory and file patterns.")
    stop("No data loaded. Check log for details.")
  }
  
  log_message(paste("Successfully loaded data from", files_loaded, "files"))
  
  
  # Set factor levels for proper ordering in plots
  log_message("Setting factor levels...")
  all_combined_lrs$relationship_type <- factor(
    all_combined_lrs$relationship_type, 
    levels = relationship_order,
    labels = relationship_labels
  )
  
  all_combined_lrs$base_loci_set <- factor(
    all_combined_lrs$base_loci_set, 
    levels = loci_set_order,
    labels = loci_set_labels
  )
  
  all_combined_lrs$population <- factor(
    all_combined_lrs$population, 
    levels = populations
  )
  
  if (nrow(all_proportions) > 0) {
    all_proportions$relationship_type <- factor(
      all_proportions$relationship_type, 
      levels = relationship_order,
      labels = relationship_labels
    )
    
    all_proportions$base_loci_set <- factor(
      all_proportions$base_loci_set, 
      levels = loci_set_order,
      labels = loci_set_labels
    )
    
    all_proportions$population <- factor(
      all_proportions$population, 
      levels = populations
    )
  }
  
  if (nrow(all_summary_stats) > 0) {
    all_summary_stats$relationship_type <- factor(
      all_summary_stats$relationship_type, 
      levels = relationship_order,
      labels = relationship_labels
    )
    
    all_summary_stats$base_loci_set <- factor(
      all_summary_stats$base_loci_set, 
      levels = loci_set_order,
      labels = loci_set_labels
    )
    
    all_summary_stats$population <- factor(
      all_summary_stats$population, 
      levels = populations
    )
  }
  
  
  # Get tallies of populations and relationship_types
  pop_rel_tallies <- all_combined_lrs %>% group_by(population, relationship_type) %>% tally()
  
  return(list(
    combined_lrs = all_combined_lrs,
    cutoffs = all_cutoffs,
    proportions = all_proportions,
    summary_stats = all_summary_stats,
    population_relationship_tallies = pop_rel_tallies
  ))
}

# Function to create boxplot of LR distributions
plot_lr_distributions <- function(combined_lrs, population_relationship_tallies) {
  log_message("Creating boxplot of LR distributions...")
  
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
  names(pop_colors) <- levels(combined_lrs$population)
  
  tryCatch({
    p <- ggplot(combined_lrs, aes(x = relationship_type, y = LR, fill = population)) +
      geom_boxplot(position = position_dodge(width = 0.85), alpha = 0.7, outlier.size = 0.5) +
      facet_wrap(~ base_loci_set, scales = "fixed") +
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
    
    ggsave(file.path(output_dir, "matched_pop_lr_boxplot.png"), plot = p, width = 14, height = 10, dpi = 300)
    log_message("LR distributions plot saved to matched_pop_lr_boxplot.png")
    
    return(p)
  }, error = function(e) {
    log_message(paste("Error creating boxplot:", e$message))
    return(NULL)
  })
}

# Function to create line plot of mean LRs
plot_mean_lr <- function(combined_lrs, summary_stats, population_relationship_tallies) {
  log_message("Creating line plot of mean LRs...")
  
  if (is.null(summary_stats) || nrow(summary_stats) == 0) {
    log_message("No summary stats data available for line plot")
    return(NULL)
  }
  
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
    group_by(relationship_type) %>%
    summarize(total_n = sum(n), .groups = 'drop')
  
  # Create facet labels with relationship types and total sample counts
  facet_labels <- setNames(
    paste0(relationship_totals$relationship_type, "\n(n=", format(relationship_totals$total_n, big.mark = ","), ")"),
    relationship_totals$relationship_type
  )
  
  # Define a color palette
  pop_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
  names(pop_colors) <- levels(summary_stats$population)
  
  tryCatch({
    p <- ggplot(summary_stats, aes(x = base_loci_set, y = mean_LR, group = population, color = population)) +
      geom_line(size = 1.2) +
      geom_point(size = 3) +
      facet_wrap(~ relationship_type, scales = "fixed", ncol = 2,
                 labeller = labeller(relationship_type = facet_labels)) +
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
    
    ggsave(file.path(output_dir, "matched_pop_mean_lr_line.png"), plot = p, width = 14, height = 10, dpi = 300)
    log_message("Mean LR plot saved to matched_pop_mean_lr_line.png")
    
    return(p)
  }, error = function(e) {
    log_message(paste("Error creating line plot:", e$message))
    return(NULL)
  })
}

# Function to create violin plots of LR distributions
plot_violin_lr <- function(combined_lrs, population_relationship_tallies) {
  log_message("Creating violin plot of LR distributions...")
  
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
  names(pop_colors) <- levels(combined_lrs$population)
  
  tryCatch({
    p <- ggplot(combined_lrs, aes(x = relationship_type, y = LR, fill = population)) +
      geom_violin(position = position_dodge(width = 0.85), alpha = 0.7, scale = "width", trim = TRUE) +
      facet_wrap(~ base_loci_set, scales = "fixed") +
      scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      scale_fill_manual(values = pop_colors, labels = legend_labels) +
      labs(
        title = "LR Distribution Densities Across Populations and Relationship Types",
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
    
    ggsave(file.path(output_dir, "matched_pop_violin_lr.png"), plot = p, width = 14, height = 10, dpi = 300)
    log_message("Violin plot saved to matched_pop_violin_lr.png")
    
    return(p)
  }, error = function(e) {
    log_message(paste("Error creating violin plot:", e$message))
    return(NULL)
  })
}

# Function to plot proportions exceeding cutoffs
plot_proportions_exceeding_cutoffs <- function(proportions_exceeding_cutoffs) {
  log_message("Creating bar plot of proportions exceeding cutoffs...")
  
  if (is.null(proportions_exceeding_cutoffs) || nrow(proportions_exceeding_cutoffs) == 0) {
    log_message("No proportions data available for bar plot")
    return(NULL)
  }
  
  # Define a color palette
  pop_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
  names(pop_colors) <- levels(proportions_exceeding_cutoffs$population)
  
  tryCatch({
    # Convert the proportions to long format
    proportions_long <- proportions_exceeding_cutoffs %>%
      pivot_longer(cols = starts_with("proportion_exceeding"),
                   names_to = "Cutoff_Type", values_to = "Proportion",
                   names_prefix = "proportion_exceeding_")
    
    proportions_long$Cutoff_Type <- factor(
      proportions_long$Cutoff_Type, 
      levels = c("fixed", "1", "0_1", "0_01"),
      labels = c("Fixed Cutoff (1.00)", "1% FPR", "0.1% FPR", "0.01% FPR")
    )
    
    p <- ggplot(proportions_long, aes(x = relationship_type, y = Proportion, fill = population)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
      facet_grid(base_loci_set ~ Cutoff_Type) +
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
    
    ggsave(file.path(output_dir, "matched_pop_proportions_exceeding_cutoffs.png"), plot = p, width = 18, height = 14, dpi = 300)
    log_message("Proportions exceeding cutoffs plot saved to matched_pop_proportions_exceeding_cutoffs.png")
    
    return(p)
  }, error = function(e) {
    log_message(paste("Error creating proportions plot:", e$message))
    return(NULL)
  })
}

# Function to create heatmap of proportions exceeding cutoffs
plot_heatmap_proportions <- function(proportions_exceeding_cutoffs) {
  log_message("Creating heatmap of proportions exceeding cutoffs...")
  
  if (is.null(proportions_exceeding_cutoffs) || nrow(proportions_exceeding_cutoffs) == 0) {
    log_message("No proportions data available for heatmap")
    return(NULL)
  }
  
  tryCatch({
    # Convert the proportions to long format
    proportions_long <- proportions_exceeding_cutoffs %>%
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
    
    p <- ggplot(proportions_long, aes(x = base_loci_set, y = relationship_type, fill = Proportion)) +
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
    
    ggsave(file.path(output_dir, "matched_pop_heatmap_proportions.png"), plot = p, width = 16, height = 12, dpi = 300)
    log_message("Heatmap of proportions exceeding cutoffs saved to matched_pop_heatmap_proportions.png")
    
    return(p)
  }, error = function(e) {
    log_message(paste("Error creating heatmap:", e$message))
    return(NULL)
  })
}

# Main execution
log_message("Starting matched population plotting process...")

# Load all data
tryCatch({
  all_data <- load_combined_data()
  
  # Create and save plots
  if (!is.null(all_data$combined_lrs) && nrow(all_data$combined_lrs) > 0) {
    lr_boxplot <- plot_lr_distributions(all_data$combined_lrs, all_data$population_relationship_tallies)
    violin_plot <- plot_violin_lr(all_data$combined_lrs, all_data$population_relationship_tallies)
  } else {
    log_message("No combined_lrs data available for plotting")
  }
  
  if (!is.null(all_data$summary_stats) && nrow(all_data$summary_stats) > 0) {
    mean_lr_plot <- plot_mean_lr(all_data$combined_lrs, all_data$summary_stats,  all_data$population_relationship_tallies)
  } else {
    log_message("No summary_stats data available for plotting")
  }
  
  if (!is.null(all_data$proportions) && nrow(all_data$proportions) > 0) {
    proportions_plot <- plot_proportions_exceeding_cutoffs(all_data$proportions)
    heatmap_plot <- plot_heatmap_proportions(all_data$proportions)
  } else {
    log_message("No proportions data available for plotting")
  }
  
  log_message("All available plots created successfully")
}, error = function(e) {
  log_message(paste("Critical error in main execution:", e$message))
  log_message("Script execution failed. Please check the logs and data files.")
  stop(paste("Critical error:", e$message))
})

log_message("Matched population plotting complete!")