#!/usr/bin/env Rscript

#### Plots of proportions exceeding cutoffs for tested parent-child and full-sibling relationships
#### Population Status: MATCHED ✓ (tested_population == population means you're using the correct allele frequencies)
#### Relationship Status: MISMATCHED ✗ (known_relationship != tested_relationship means you're testing the wrong relationship hypothesis)


# Load Required Libraries
suppressMessages(suppressWarnings({
  library(tidyverse)  # Includes ggplot2, dplyr, tidyr, etc.
  library(data.table) # Required for efficient fread()
  library(scales)     # For number formatting
  library(shades)     # For nice colors
}))

# Source required modules
source("code/module9_combinedLR_stats_functions.R")

# Helper function for logging
log_message <- function(message) {
  cat(paste0("[", Sys.time(), "] ", message, "\n"))
}

# --- Argument Parsing and Setup ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript code/plots_mismatched.R <input_dir> [output_dir]")
}

input_dir <- args[1]
log_message(paste("Input directory:", input_dir))

if (length(args) >= 2) {
  output_subdir <- args[2]
  output_dir <- file.path("output", output_subdir)
} else {
  timestamp <- format(Sys.time(), "%Y%m%d")
  output_dir <- file.path(input_dir, paste0("mismatched_pop_plots_", timestamp)) 
}
log_message(paste("Output directory:", output_dir))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)



today <- format(Sys.time(), format = "%Y-%m-%d")


# Define factors, colors, and labels
populations <- c("AfAm", "Cauc", "Hispanic", "Asian", "all")

relationship_order <- c("parent_child", "full_siblings", "half_siblings", 
                        "cousins", "second_cousins", "unrelated")
relationship_labels <- c("Parent-Child", "Full Siblings", "Half Siblings", 
                         "Cousins", "Second Cousins", "Unrelated")
names(relationship_labels) <- relationship_order

loci_set_order <- c("core_13", "identifiler_15", "expanded_20", "supplementary", "autosomal_29")
loci_set_labels <- c("Core 13", "Identifiler 15", "Expanded 20", "Supplementary", "Autosomal 29")
names(loci_set_labels) <- loci_set_order

# Define color palette for populations (frequency sources)
true_pop_colors <- c(
  "AfAm" = "#E41A1C",     # Red
  "Asian" = "#377EB8",    # Blue
  "Cauc" = "#4DAF4A",     # Green
  "Hispanic" = "#984EA3",  # Purple
  "all" = "#FF7F00"       # Orange
)

unnamed_colors <- unname(true_pop_colors)
light_colors_unnamed <- lightness(unnamed_colors, scalefac(1.2))
light_pop_colors <- light_colors_unnamed
names(light_pop_colors) <- names(true_pop_colors)


######## Read in combined_LR files ########
# combined_LR_all.rds made in analyze_lr_outputs.R
all_combined_file <- file.path(output_dir, "combined_LR_all.rds")
all_combined <- readRDS(all_combined_file)

######## Make sure LR is numeric ######## 
all_combined <- all_combined %>% mutate(across(c(combined_LR), as.numeric))
all_combined <- all_combined %>% mutate(
  across(c(known_relationship, tested_relationship), 
         ~factor(., levels = relationship_order)),
  loci_set = factor(loci_set, levels = loci_set_order)
)

cutoffs <- calculate_cutoffs(all_combined, c(1, 0.1, 0.01))
cutoffs_file <- file.path(output_dir, "combined_LR_cutoffs.csv")
fwrite(cutoffs, cutoffs_file)

proportions_exceeding_cutoffs <- calculate_proportions_exceeding_cutoffs(all_combined, cutoffs)
proportions_file <- file.path(output_dir, "combined_LR_exceeding_cutoffs.csv")
fwrite(proportions_exceeding_cutoffs, proportions_file)


######## Bar chart plots ######## 
# Open PDF device to save all plots to one PDF
pdf(file.path(output_dir, "population_mismatch_proportions_analysis.pdf"), width = 14, height = 10)

# For relationship mismatch analysis (mimicking the single image)
relationships_to_test <- c("parent_child", "full_siblings")

for (rel in relationships_to_test) {
  # Filter to include ALL relationships (both correct and incorrect)
  rel_data <- proportions_exceeding_cutoffs %>% 
    filter(tested_relationship == rel,           # Testing this relationship
           tested_population == population)      # Using correct population frequencies
  
  # Get counts
  counts_summary <- rel_data %>%
    group_by(population, known_relationship, loci_set) %>%
    summarize(n_pairs = first(n_related), .groups = 'drop')
  
  cat("For", rel, " hypothesis each population x known_relationship x loci set has:", unique(counts_summary$n_pairs), "pairs\n")
  
  # Create long format data for plotting
  proportions_long <- rel_data %>%
    pivot_longer(cols = starts_with("proportion_exceeding"),
                 names_to = "Cutoff_Type", values_to = "Proportion",
                 names_prefix = "proportion_exceeding_")
  
  proportions_long$Cutoff_Type <- factor(proportions_long$Cutoff_Type, 
                                         levels = c("0_01", "0_1","1", "fixed"),
                                         labels = c( "0.01% FPR", "0.1% FPR", "1% FPR","Fixed (1.00)"))
  
  proportions_long <- proportions_long %>%
    mutate(
      known_relationship = recode_factor(known_relationship,
                                         "parent_child" = "Parent-Child",
                                         "full_siblings" = "Full Siblings", 
                                         "half_siblings" = "Half Siblings",
                                         "cousins" = "Cousins",
                                         "second_cousins" = "Second Cousins",
                                         "unrelated" = "Unrelated"),
      loci_set = recode_factor(loci_set,
                               "core_13" = "Core 13",
                               "identifiler_15" = "Identifiler 15",
                               "expanded_20" = "Expanded 20", 
                               "supplementary" = "Supplementary",
                               "autosomal_29" = "Autosomal 29")
    )
  
  
  # Plot with known_relationship on x-axis (this creates the variation!)
  p <- ggplot(proportions_long, aes(x = known_relationship, y = Proportion, fill = population)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    facet_grid(Cutoff_Type ~ loci_set, scales = "fixed", labeller = labeller(loci_set = loci_set_labels)) +
    scale_fill_manual(values = true_pop_colors, name = "True Population") +
    labs(
      title = paste("Proportions exceeding likelihood cutoffs:", rel, "Hypothesis"),
      subtitle = paste("Using population-matched allele frequencies"),
      x = "True Relationship",
      y = "Proportion Exceeding Cutoff"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 14, color = "darkblue"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.title = element_text(size = 12),
      strip.text = element_text(size = 9),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )
  
  print(p)
  
  # Save individual PNG
  filename <- paste0("relationship_analysis_", rel, "_hypot_", today, ".png")
  ggsave(file.path(output_dir, filename), plot = p, width = 14, height = 10, bg = "white", dpi = 300)
}

dev.off()