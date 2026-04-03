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

    log10_LR = suppressWarnings(log10(combined_LR)),
    
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
    median_log10_LR = median(log10_LR, na.rm = TRUE),
    mean_log10_LR   = mean(log10_LR[is.finite(log10_LR)], na.rm = TRUE),
    se_log10_LR     = sd(log10_LR[is.finite(log10_LR)], na.rm = TRUE) / sqrt(sum(combined_LR > 0)),
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
# SUPPLEMENTARY FIGURE: Population Mismatch Percent Change Heatmap
# Shows % change in LR when using wrong vs. correct population frequencies,
# computed per pair as log10(LR_wrong / LR_correct), then summarized as
# median across pairs. Diagonal is computed empirically (should be ~0%).
# ==============================================================================

log_message("Generating Supplementary Figure: Population Mismatch Heatmap...")

# Compute median log10(LR) per population x tested_population x relationship x loci cell
# Zeros are included: log10(0) = -Inf, median handles correctly since <50% per group are zero
correct_lrs <- all_combined %>%
  filter(pop_match == TRUE, rel_match == TRUE,
         known_relationship %in% c("Parent-Child", "Full Siblings")) %>%
  select(batch_id, pair_id, population, known_relationship,
         loci_set, tested_relationship,
         correct_log10_LR = log10_LR)

wrong_lrs <- all_combined %>%
  filter(pop_match == FALSE, rel_match == TRUE,
         known_relationship %in% c("Parent-Child", "Full Siblings")) %>%
  select(batch_id, pair_id, population, known_relationship,
         loci_set, tested_relationship, tested_population,
         wrong_log10_LR = log10_LR)

figS_data <- wrong_lrs %>%
  left_join(correct_lrs,
            by = c("batch_id", "pair_id", "population",
                   "known_relationship", "loci_set", "tested_relationship")) %>%
  mutate(log10_ratio = wrong_log10_LR - correct_log10_LR) %>%
  group_by(population, tested_population, known_relationship, loci_set) %>%
  summarize(
    median_log10_ratio = median(log10_ratio, na.rm = TRUE),
    n_pairs            = n(),
    .groups = "drop"
  ) %>%
  mutate(
    tested_population = factor(tested_population,
                               levels = c("AfAm", "Asian", "Cauc", "Hispanic", "all")),
    population        = factor(population,
                               levels = c("AfAm", "Asian", "Cauc", "Hispanic", "all")),
    text_color = if_else(median_log10_ratio > 1.5, "white", "black"),
    is_diagonal = FALSE  # all rows here are off-diagonal by construction; diagonal added below
  )

diagonal_rows <- all_combined %>%
  filter(pop_match == TRUE, rel_match == TRUE,
         known_relationship %in% c("Parent-Child", "Full Siblings")) %>%
  select(batch_id, pair_id, population, known_relationship,
         loci_set, tested_relationship, tested_population,
         log10_LR) %>%
  left_join(
    all_combined %>%
      filter(pop_match == TRUE, rel_match == TRUE,
             known_relationship %in% c("Parent-Child", "Full Siblings")) %>%
      select(batch_id, pair_id, population, known_relationship,
             loci_set, tested_relationship,
             correct_log10_LR = log10_LR),
    by = c("batch_id", "pair_id", "population",
           "known_relationship", "loci_set", "tested_relationship")
  ) %>%
  mutate(log10_ratio = log10_LR - correct_log10_LR) %>%
  group_by(population, tested_population, known_relationship, loci_set) %>%
  summarize(
    median_log10_ratio = median(log10_ratio, na.rm = TRUE),
    n_pairs            = n(),
    .groups = "drop"
  ) %>%
  mutate(
    tested_population = factor(tested_population,
                               levels = c("AfAm", "Asian", "Cauc", "Hispanic", "all")),
    population        = factor(population,
                               levels = c("AfAm", "Asian", "Cauc", "Hispanic", "all")),
    text_color  = "black",
    is_diagonal = TRUE
    )

figS_data <- bind_rows(figS_data, diagonal_rows)

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
