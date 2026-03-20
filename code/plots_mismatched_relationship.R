#!/usr/bin/env Rscript
################################################################################
# Publication-Ready Mismatched Relationship Analysis Plots
#
# Purpose: Generate figure examining LR distributions when pairs are tested
#          under correct vs. incorrect relationship hypotheses
#          (using correct population frequencies throughout)
#
#   Main Figure: LR distributions for all true relationship types tested
#                under parent-child and full-sibling hypotheses, across
#                Core 13, Expanded 20, and Autosomal 29 loci panels
#                (gold shading highlights correct hypothesis matches)
#
# Note: This analysis is partially redundant with the FPR cutoff analysis
#       (plots_proportion_exceeding_cutoffs.R), which summarises the same
#       patterns as binary classification rates. This script is retained
#       because it shows the full continuous LR distributions rather than
#       a single threshold-based summary statistic.
#
# Date: 2026-03-06
################################################################################

# ==============================================================================
# SETUP AND DEPENDENCIES
# ==============================================================================

suppressMessages({
  library(tidyverse)    # Data manipulation and plotting
  library(data.table)   # Fast data reading
  library(scales)       # Scale functions for plotting
  library(patchwork)    # Combining plots
})

# Helper function for timestamped logging
log_message <- function(message) {
  cat(sprintf("[%s] %s\n", Sys.time(), message))
}

# ==============================================================================
# ARGUMENT PARSING AND DIRECTORY SETUP
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript plots_mismatched_relationship.R <input_dir> [output_dir]")
}

input_dir <- args[1]
output_dir <- if (length(args) >= 2) args[2] else file.path(input_dir, "publication_figures")

log_message(sprintf("Input directory: %s", input_dir))
log_message(sprintf("Output directory: %s", output_dir))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# CONSTANTS AND STYLING
# ==============================================================================

# Define relationship ordering and labels
RELATIONSHIP_ORDER <- c("parent_child", "full_siblings", "half_siblings",
                        "cousins", "second_cousins", "unrelated")
RELATIONSHIP_LABELS <- c("Parent-Child", "Full Siblings", "Half Siblings",
                         "Cousins", "Second Cousins", "Unrelated")

# Define loci set ordering and labels
LOCI_ORDER <- c("core_13", "identifiler_15", "expanded_20",
                "supplementary", "autosomal_29")
LOCI_LABELS <- c("Core 13", "Identifiler 15", "Expanded 20",
                 "Supplementary", "Autosomal 29")

# Population colors (consistent with existing analyses)
POP_COLORS <- c(
  "AfAm"     = "#E41A1C",   # Red
  "Asian"    = "#377EB8",   # Blue
  "Cauc"     = "#4DAF4A",   # Green
  "Hispanic" = "#984EA3",   # Purple
  "all"      = "#FF7F00"    # Orange
)

# Publication theme
theme_publication <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      # Text
      plot.title    = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
      plot.subtitle = element_text(size = rel(1), hjust = 0.5, color = "gray30"),
      axis.title    = element_text(face = "bold", size = rel(1)),
      axis.text     = element_text(size = rel(0.9)),

      # Facets
      strip.background = element_rect(fill = "gray90", color = "gray60"),
      strip.text       = element_text(face = "bold", size = rel(0.9)),

      # Legend
      legend.title    = element_text(face = "bold", size = rel(1)),
      legend.text     = element_text(size = rel(0.9)),
      legend.position = "bottom",
      legend.box      = "horizontal",

      # Grid
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),

      # Margins
      plot.margin = margin(10, 10, 10, 10)
    )
}

# ==============================================================================
# DATA LOADING AND PREPROCESSING
# ==============================================================================

log_message("Loading combined LR data...")

all_combined_file <- file.path(input_dir, "combined_LR_all.rds")

if (!file.exists(all_combined_file)) {
  stop(sprintf("Data file not found: %s", all_combined_file))
}

all_combined <- readRDS(all_combined_file)

log_message(sprintf("Loaded %s rows", format(nrow(all_combined), big.mark = ",")))

# Preprocess data
all_combined <- all_combined %>%
  mutate(
    combined_LR = as.numeric(combined_LR),

    known_relationship = factor(known_relationship,
                                levels = RELATIONSHIP_ORDER,
                                labels = RELATIONSHIP_LABELS),
    tested_relationship = factor(tested_relationship,
                                 levels = RELATIONSHIP_ORDER,
                                 labels = RELATIONSHIP_LABELS),
    loci_set   = factor(loci_set,
                        levels = LOCI_ORDER,
                        labels = LOCI_LABELS),
    population = factor(population,
                        levels = c("AfAm", "Asian", "Cauc", "Hispanic", "all")),

    # Population match indicator
    pop_match       = tested_population == population,
    pop_match_label = if_else(pop_match, "Matched", "Mismatched"),

    # Relationship match indicator
    rel_match = as.character(known_relationship) == as.character(tested_relationship)
  ) %>%
  filter(!is.na(combined_LR), is.finite(combined_LR), combined_LR > 0)

log_message(sprintf("After preprocessing: %s rows", format(nrow(all_combined), big.mark = ",")))

# ==============================================================================
# MAIN FIGURE: RELATIONSHIP DISCRIMINATION
# Question: Can the method identify the correct relationship type?
# Approach: Show LR distributions when testing pairs under different hypotheses
#           Using CORRECT population frequencies only (pop_match == TRUE)
#           Gold shading highlights panels where tested hypothesis == true relationship
# ==============================================================================

log_message("Generating Main Figure: Relationship Discrimination...")

# Filter to matched populations and subset of loci panels for clarity
fig_rel_data <- all_combined %>%
  filter(
    pop_match == TRUE,
    tested_relationship %in% c("Parent-Child", "Full Siblings"),
    loci_set %in% c("Core 13", "Expanded 20", "Autosomal 29")
  )

fig_rel <- ggplot(fig_rel_data,
                  aes(x = tested_relationship, y = combined_LR,
                      fill = population)) +

  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3,
               position = position_dodge(width = 0.85), alpha = 0.8) +

  # Gold shading on correct hypothesis panels (where tested == known)
  geom_rect(
    data = fig_rel_data %>%
      filter(rel_match) %>%
      distinct(known_relationship, tested_relationship, loci_set),
    aes(xmin = as.numeric(tested_relationship) - 0.45,
        xmax = as.numeric(tested_relationship) + 0.45,
        ymin = 1, ymax = 1e40),
    fill = "gold", alpha = 0.15, inherit.aes = FALSE
  ) +

  facet_grid(known_relationship ~ loci_set, scales = "fixed") +

  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_fill_manual(values = POP_COLORS, name = "True Population") +

  labs(
    x = "Tested Relationship Hypothesis",
    y = "Combined Likelihood Ratio (LR)"
  ) +

  theme_publication() +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

fig_rel_filename <- "mismatched_relationship_discrimination.png"
ggsave(
  filename = file.path(output_dir, fig_rel_filename),
  plot     = fig_rel,
  width    = 12,
  height   = 10,
  dpi      = 300,
  bg       = "white"
)

log_message(sprintf("Main figure saved: %s", fig_rel_filename))

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\nOUTPUT FILES:\n")
cat(sprintf("   - %s (main figure)\n", fig_rel_filename))
cat("\nAll files saved to:", output_dir, "\n")
cat("================================================================================\n")

log_message("Analysis complete!")
