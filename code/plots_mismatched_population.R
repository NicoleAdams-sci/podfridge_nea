#!/usr/bin/env Rscript
################################################################################
# Publication-Ready Mismatched Population Analysis Plots
#
# Purpose: Generate figures examining the effect of using wrong population
#          allele frequencies on LR performance (population robustness)
#
#   Main Figure:        LR comparison for matched vs. mismatched population
#                       frequencies across loci panels (true positives only)
#   Supplementary Figure: heatmap of median log10(LR) inflation due to
#                         population frequency mismatch, across all loci panels
#   Summary Table:      Detailed population comparison statistics
#
# Input:  <input_dir>/mismatched_pop_robustness.csv
#         <input_dir>/mismatched_pop_heatmap.csv
#         (produced by prepare_combined_lr_intermediates.R)
#
# Usage:
#   Rscript code/plots_mismatched_population.R <input_dir> [output_dir]
#
#   input_dir   Full path to data directory (e.g., output/lr_analysis_20260130)
#   output_dir  Where to write figures (default: <input_dir>/plots_mismatched_population)
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
  stop("Usage: Rscript code/plots_mismatched_population.R <input_dir> [output_dir]")
}

input_dir  <- args[1]
output_dir <- if (length(args) >= 2) args[2] else file.path(input_dir, "plots_mismatched_population")

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

# Classification colors
MATCH_COLORS <- c(
  "Matched"    = "#2166ac",   # Blue
  "Mismatched" = "#d6604d"    # Red
)

# Publication theme
theme_publication <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(
      # Text
      plot.title    = element_text(size = rel(1.2), hjust = 0.5),
      plot.subtitle = element_text(size = rel(1), hjust = 0.5, color = "gray30"),
      axis.title    = element_text( size = rel(1)),
      axis.text     = element_text(size = rel(0.9)),

      # Facets
      strip.background = element_rect(fill = "gray90", color = "gray60"),
      strip.text       = element_text(size = rel(0.8)),

      # Legend
      legend.title    = element_text(size = rel(1)),
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
# DATA LOADING
# ==============================================================================

log_message("Loading intermediate files...")

robustness_file <- file.path(input_dir, "mismatched_pop_robustness.csv")
heatmap_file    <- file.path(input_dir, "mismatched_pop_heatmap.csv")

for (f in c(robustness_file, heatmap_file)) {
  if (!file.exists(f)) {
    stop(sprintf(
      "Intermediate file not found: %s\nRun prepare_combined_lr_intermediates.R first.", f
    ))
  }
}

fig_pop_data <- fread(robustness_file) %>%
  mutate(
    known_relationship = factor(known_relationship,
                                levels = RELATIONSHIP_ORDER,
                                labels = RELATIONSHIP_LABELS),
    loci_set           = factor(loci_set,
                                levels = LOCI_ORDER,
                                labels = LOCI_LABELS),
    population         = factor(population,
                                levels = c("AfAm", "Asian", "Cauc", "Hispanic", "all")),
    pop_match_label    = factor(pop_match_label, levels = c("Matched", "Mismatched"))
  )

figS_data <- fread(heatmap_file) %>%
  mutate(
    known_relationship = factor(known_relationship,
                                levels = RELATIONSHIP_ORDER,
                                labels = RELATIONSHIP_LABELS),
    loci_set           = factor(loci_set,
                                levels = LOCI_ORDER,
                                labels = LOCI_LABELS),
    tested_population  = factor(tested_population,
                                levels = c("AfAm", "Asian", "Cauc", "Hispanic", "all")),
    population         = factor(population,
                                levels = c("AfAm", "Asian", "Cauc", "Hispanic", "all")),
    text_color         = if_else(median_log10_ratio > 1.5, "white", "black"),
    # diagonal rows have median_log10_ratio = 0; force text_color to black
    text_color         = if_else(is_diagonal, "black", text_color)
  )

log_message(sprintf("Loaded %d rows (robustness), %d rows (heatmap).",
                    nrow(fig_pop_data), nrow(figS_data)))

fig_pop <- ggplot(fig_pop_data,
                  aes(x = n_loci, y = median_log10_LR,
                      color    = pop_match_label,
                      linetype = pop_match_label,
                      group    = interaction(population, tested_population,
                                             pop_match_label))) +

  geom_line(alpha = 0.6, linewidth = 0.8) +
  geom_point(size = 2, alpha = 0.8) +

  facet_grid(population ~ known_relationship) +

  scale_x_continuous(breaks = c(13, 15, 20, 23, 29)) +
  scale_color_manual(values = MATCH_COLORS,
                     name = "Population\nFrequencies") +
  scale_linetype_manual(values = c("Matched" = "solid", "Mismatched" = "dashed"),
                        name = "Population\nFrequencies") +

  labs(
    x = "Number of Loci",
    y = "Median log\u2081\u2080(Combined LR)"
  ) +

  theme_publication() +
  theme(
    legend.position = "right",
    plot.caption    = element_text(hjust = 1, face = "italic", color = "gray50")
  )

fig_pop_filename <- "mismatched_population_robustness.png"
ggsave(
  filename = file.path(output_dir, fig_pop_filename),
  plot     = fig_pop,
  width    = 10,
  height   = 8,
  dpi      = 300,
  bg       = "white"
)

log_message(sprintf("Main figure saved: %s", fig_pop_filename))

# ==============================================================================
# SUPPLEMENTARY FIGURE: Population Mismatch Percent Change Heatmap
# Shows % change in LR when using wrong vs. correct population frequencies,
# computed per pair as log10(LR_wrong / LR_correct), then summarized as
# median across pairs. Diagonal is computed empirically (should be ~0%).
# ==============================================================================

log_message("Generating Supplementary Figure: Population Mismatch Heatmap...")

figS_pop <- ggplot(figS_data,
                   aes(x = tested_population, y = population,
                       fill = median_log10_ratio)) +

  geom_tile(color = "white", linewidth = 1) +

  geom_text(aes(label = sprintf("%+.2f", median_log10_ratio)),
                color = figS_data$text_color, size = 3.8) +

  # Highlight diagonal (delta = 0)
  geom_tile(data = filter(figS_data, is_diagonal),
            color = "gold", linewidth = 2, fill = NA) +

  # Rows = loci panels, columns = relationship type
  facet_grid(loci_set ~ known_relationship) +

  # scale_fill_gradientn(
  #   colors = c("white", "#c6dbef", "#6baed6", "#2166ac", "#08306b"),
  #   trans  = "log1p",   # log1p handles 0 gracefully
  #   name   = "Median\nlog\u2081\u2080(LR_wrong /\nLR_correct)"
  # ) +
  
  scale_fill_gradientn(
    colors = c("white", "#c6dbef", "#6baed6", "#2166ac", "#08306b"),
    limits = c(0, 3),
    oob    = scales::squish,   # values above 3 get capped to darkest color
    name = "Median\nlog\u2081\u2080(LR_wrong /\nLR_correct)\n(capped at 3)"
  ) +

  labs(
    title    = "Population Frequency Mismatch: LR Inflation Relative to Matched Frequencies",
    subtitle = "Each cell = median log\u2081\u2080(LR_wrong / LR_correct) across all pairs\n+1 = wrong frequencies give 10\u00d7 higher LR | 0 = no effect (gold border)",
    x        = "Tested Population Frequencies",
    y        = "True Population"
  ) +

  theme_publication() +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    panel.grid      = element_blank(),
    legend.position = "right"
  ) +

  coord_equal()

figS_filename <- "mismatched_population_supp_heatmap.png"
ggsave(
  filename = file.path(output_dir, figS_filename),
  plot     = figS_pop,
  width    = 12,
  height   = 20,
  dpi      = 300,
  bg       = "white"
)

log_message(sprintf("Supplementary figure saved: %s", figS_filename))

# ==============================================================================
# SUMMARY STATISTICS TABLE
# ==============================================================================

log_message("Generating summary table...")

detailed_stats <- fig_pop_data %>%
  select(population, tested_population, known_relationship, loci_set,
         pop_match_label, median_log10_LR, mean_log10_LR, se_log10_LR,
         n_pairs, n_zero) %>%
  arrange(population, known_relationship, loci_set, tested_population)

detailed_filename <- "mismatched_population_comparison_detailed.csv"
write_csv(detailed_stats, file.path(output_dir, detailed_filename))

log_message(sprintf("Summary table saved: %s", detailed_filename))

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\nOUTPUT FILES:\n")
cat(sprintf("   - %s (main figure)\n",        fig_pop_filename))
cat(sprintf("   - %s (supplementary figure)\n", figS_filename))
cat(sprintf("   - %s (summary table)\n",       detailed_filename))
cat("\nAll files saved to:", output_dir, "\n")
cat("================================================================================\n")

log_message("Analysis complete!")
