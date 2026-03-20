#!/usr/bin/env Rscript
################################################################################
# Publication-Ready Mismatched Population Analysis Plots
#
# Purpose: Generate figures examining the effect of using wrong population
#          allele frequencies on LR performance (population robustness)
#
#   Main Figure:        LR comparison for matched vs. mismatched population
#                       frequencies across loci panels (true positives only)
#   Supplementary Figure: Delta heatmap of median log10(LR) inflation due to
#                         population frequency mismatch, across all loci panels
#   Summary Table:      Detailed population comparison statistics
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
  stop("Usage: Rscript plots_mismatched_population.R <input_dir> [output_dir]")
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

# Classification colors
MATCH_COLORS <- c(
  "Matched"    = "#2166ac",   # Blue
  "Mismatched" = "#d6604d"    # Red
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

# ==============================================================================
# DIAGNOSTIC: CHECK FOR ZERO / NON-FINITE LRs (PRE-FILTER)
# Run before filtering so we capture what is being dropped.
# Results are printed to the SLURM log for post-hoc inspection.
# ==============================================================================

log_message("Running LR diagnostic checks (pre-filter)...")

all_combined <- all_combined %>%
  mutate(combined_LR = as.numeric(combined_LR))

n_total_raw  <- nrow(all_combined)
n_zero_raw   <- sum(all_combined$combined_LR == 0,   na.rm = TRUE)
n_na_raw     <- sum(is.na(all_combined$combined_LR))
n_neginf_raw <- sum(!is.finite(log10(suppressWarnings(
                  all_combined$combined_LR))), na.rm = TRUE)

cat(sprintf("\n--- LR Diagnostic (pre-filter, n = %s) ---\n",
            format(n_total_raw, big.mark = ",")))
cat(sprintf("  Zero LRs:              %d (%.4f%%)\n",
            n_zero_raw,   100 * n_zero_raw   / n_total_raw))
cat(sprintf("  NA LRs:                %d (%.4f%%)\n",
            n_na_raw,     100 * n_na_raw     / n_total_raw))
cat(sprintf("  Non-finite log10(LR):  %d (%.4f%%)\n",
            n_neginf_raw, 100 * n_neginf_raw / n_total_raw))

# Exclusion rate by population match — key check for whether mismatch
# drives more exclusions than matched frequencies
excl_by_match <- all_combined %>%
  mutate(
    is_zero_or_na = is.na(combined_LR) | combined_LR == 0,
    pop_match     = tested_population == population
  ) %>%
  group_by(pop_match) %>%
  summarize(
    n_total      = n(),
    n_excluded   = sum(is_zero_or_na),
    pct_excluded = 100 * sum(is_zero_or_na) / n(),
    .groups = "drop"
  )

cat("\n--- Exclusion rate by population match (pre-filter) ---\n")
print(as.data.frame(excl_by_match))
cat("\n")

# Break down zero LRs by known vs. tested relationship combination
# Identifies whether zeros cluster in specific relationship mismatches
# (e.g., unrelated pairs tested as parent-child where k0 = 0)
zero_by_relationship <- all_combined %>%
  filter(combined_LR == 0) %>%
  count(known_relationship, tested_relationship, name = "n_zeros") %>%
  arrange(desc(n_zeros)) %>%
  mutate(pct_of_all_zeros = 100 * n_zeros / n_zero_raw)

cat("--- Zero LR breakdown by known vs. tested relationship (pre-filter) ---\n")
print(as.data.frame(zero_by_relationship))
cat("\n")

# ==============================================================================
# PREPROCESS DATA
# ==============================================================================

all_combined <- all_combined %>%
  mutate(
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
  # Keep zeros — legitimate exclusions under parent-child test (k0=0)
  # Median is robust to these as long as <50% per group are zero
  filter(!is.na(combined_LR))

log_message(sprintf("After preprocessing: %s rows", format(nrow(all_combined), big.mark = ",")))

# Post-filter diagnostic
n_total_clean <- nrow(all_combined)
n_zero_clean  <- sum(all_combined$combined_LR == 0, na.rm = TRUE)

cat(sprintf("\n--- LR Diagnostic (post-filter, n = %s) ---\n",
            format(n_total_clean, big.mark = ",")))
cat(sprintf("  Zero LRs retained:     %d (%.4f%%)\n\n",
            n_zero_clean, 100 * n_zero_clean / n_total_clean))

# ==============================================================================
# MAIN FIGURE: POPULATION ROBUSTNESS
# Question: How sensitive is the method to using wrong population frequencies?
# Approach: Compare LRs for TRUE POSITIVE pairs when using matched vs.
#           mismatched population frequency tables
# ==============================================================================

log_message("Generating Main Figure: Population Robustness...")

# Filter to true positives only (relationship correctly tested)
fig_pop_data <- all_combined %>%
  filter(
    rel_match == TRUE,
    known_relationship %in% c("Parent-Child", "Full Siblings")
  ) %>%
  group_by(population, known_relationship, loci_set,
           pop_match_label, tested_population) %>%
  summarize(
    # Median is robust to zero LRs (legitimate parent-child exclusions)
    # log10(0) = -Inf; median handles this correctly since <50% per group are zero
    median_log10_LR = median(suppressWarnings(log10(combined_LR)), na.rm = TRUE),
    # Mean and SE computed on non-zero pairs only — retained for summary CSV
    mean_log10_LR   = mean(suppressWarnings(log10(combined_LR[combined_LR > 0])),
                           na.rm = TRUE),
    se_log10_LR     = sd(suppressWarnings(log10(combined_LR[combined_LR > 0])),
                         na.rm = TRUE) / sqrt(sum(combined_LR > 0)),
    n_pairs         = n(),
    n_zero          = sum(combined_LR == 0),
    .groups = "drop"
  ) %>%
  mutate(
    n_loci = case_when(
      loci_set == "Core 13"        ~ 13,
      loci_set == "Identifiler 15" ~ 15,
      loci_set == "Expanded 20"    ~ 20,
      loci_set == "Supplementary"  ~ 23,
      loci_set == "Autosomal 29"   ~ 29
    ),
    pop_match_label = factor(pop_match_label, levels = c("Matched", "Mismatched"))
  )

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
# SUPPLEMENTARY FIGURE: Population Mismatch Delta Heatmap
# Shows delta in median log10(LR) between mismatched and matched population
# frequencies for all pairwise population combinations, across all loci panels
# (rows) and relationship types (columns).
# Diagonal = 0 by definition (matched - matched).
# Positive delta = mismatch inflates LR relative to correct frequencies.
# Includes "all" pooled population as reference.
# ==============================================================================

log_message("Generating Supplementary Figure: Population Mismatch Delta Heatmap...")

# Compute median log10(LR) per population x tested_population x relationship x loci cell
# Zeros are included: log10(0) = -Inf, median handles correctly since <50% per group are zero
figS_raw <- all_combined %>%
  filter(
    rel_match == TRUE,
    known_relationship %in% c("Parent-Child", "Full Siblings")
  ) %>%
  group_by(population, tested_population, known_relationship, loci_set) %>%
  summarize(
    median_log10_LR = median(suppressWarnings(log10(combined_LR)), na.rm = TRUE),
    n_pairs         = n(),
    .groups = "drop"
  ) %>%
  # Ensure "all" appears last on both axes
  mutate(
    tested_population = factor(tested_population,
                               levels = c("AfAm", "Asian", "Cauc", "Hispanic", "all")),
    population        = factor(population,
                               levels = c("AfAm", "Asian", "Cauc", "Hispanic", "all"))
  )

# Extract matched (diagonal) median for each population x relationship x loci_set
matched_vals <- figS_raw %>%
  filter(population == tested_population) %>%
  select(population, known_relationship, loci_set,
         matched_median = median_log10_LR)

# Compute delta: mismatched median - matched median
# Diagonal will be 0 by definition
figS_data <- figS_raw %>%
  left_join(matched_vals,
            by = c("population", "known_relationship", "loci_set")) %>%
  mutate(
    delta           = median_log10_LR - matched_median,
    is_diagonal     = population == tested_population
  )

figS_pop <- ggplot(figS_data,
                   aes(x = tested_population, y = population,
                       fill = delta)) +

  geom_tile(color = "white", linewidth = 1) +

  geom_text(aes(label = sprintf("%+.1f", delta)),
            color = "black", size = 2.8, fontface = "bold") +

  # Highlight diagonal (delta = 0)
  geom_tile(data = filter(figS_data, is_diagonal),
            color = "gold", linewidth = 2, fill = NA) +

  # Rows = loci panels, columns = relationship type
  facet_grid(loci_set ~ known_relationship) +

  scale_fill_gradient(
    low    = "white",
    high   = "#2166ac",   # Blue = higher inflation
    name   = "\u0394 Median\nlog\u2081\u2080(LR)\n(Mismatch \u2212 Matched)",
    breaks = pretty_breaks(n = 5)
  ) +

  labs(
    title    = "Population Frequency Mismatch: LR Inflation Relative to Matched Frequencies",
    subtitle = "Delta = median log\u2081\u2080(LR) using wrong frequencies \u2212 median log\u2081\u2080(LR) using correct frequencies\nGold borders = correct population match (delta = 0)",
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

figS_filename <- "mismatched_population_supp_delta_heatmap.png"
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
