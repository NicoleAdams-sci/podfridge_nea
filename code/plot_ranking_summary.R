#!/usr/bin/env Rscript
################################################################################
# Focal Ranking Test - Publication-Ready Summary Plots
#
# Purpose: Summarize the focal ranking test - for simulated true relatives
#          (parent-child, full siblings, half siblings) ranked against a
#          pooled unrelated database, how often does the true relative land
#          in the top N candidates, and how far does rank degrade when the
#          pair is tested under the WRONG relationship hypothesis?
#
#   Figure 1: Proportion of true relatives recovered in top 10/50/100/200,
#             faceted by true relationship (rows) x loci set (columns),
#             colored by tested hypothesis - gold border highlights bars
#             where the tested hypothesis matches the true relationship
#   Figure 2: Distribution of true relative rank (log10), same faceting,
#             gold shading highlights matched-hypothesis panels
#   Figure 3 (optional): Tied-group size - how many candidates shared the
#             true relative's exact combined_LR (large values mean rank is
#             an average over an indistinguishable block, not a real
#             position - common when combined_LR == 0 under a mismatched
#             hypothesis). Only produced if the input file has this column.
#   Summary table: descriptive statistics for the results section
#
# Input:  combined ranking_outcomes CSV (see combine_ranking_outcomes.sh)
#
# Usage:
#   Rscript code/plot_ranking_summary.R <combined_outcomes.csv> [output_dir]
#
# Date: 2026-07-09
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

log_message <- function(message) {
  cat(sprintf("[%s] %s\n", Sys.time(), message))
}

# ==============================================================================
# ARGUMENT PARSING AND DIRECTORY SETUP
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript plot_ranking_summary.R <combined_outcomes.csv> [output_dir]")
}

input_file <- args[1]
output_dir <- if (length(args) >= 2) args[2] else file.path(dirname(input_file), "publication_figures")

log_message(sprintf("Input file: %s", input_file))
log_message(sprintf("Output directory: %s", output_dir))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# CONSTANTS AND STYLING
# ==============================================================================

# Relationship ordering/labels (consistent with plots_matched_publication.R,
# plots_mismatched_relationship.R, plots_mismatched_population.R)
RELATIONSHIP_ORDER <- c("parent_child", "full_siblings", "half_siblings",
                        "cousins", "second_cousins", "unrelated")
RELATIONSHIP_LABELS <- c("Parent-Child", "Full Siblings", "Half Siblings",
                         "Cousins", "Second Cousins", "Unrelated")

# Loci set ordering/labels
LOCI_ORDER <- c("core_13", "identifiler_15", "expanded_20",
                "supplementary", "autosomal_29")
LOCI_LABELS <- c("Core 13", "Identifiler 15", "Expanded 20",
                 "Supplementary", "Autosomal 29")

# Relationship colors (Okabe-Ito colorblind-safe palette, shared across scripts)
relationship_colors <- c(
  "Parent-Child"   = "#D55E00",   # Vermillion
  "Full Siblings"  = "#E69F00",   # Orange
  "Half Siblings"  = "#56B4E9",   # Sky blue
  "Cousins"        = "#009E73",   # Bluish green
  "Second Cousins" = "#CC79A7",   # Reddish purple
  "Unrelated"      = "#999999"    # Gray
)

# Publication theme (matches plots_mismatched_population.R)
theme_publication <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title    = element_text(size = rel(1.2), hjust = 0.5),
      plot.subtitle = element_text(size = rel(0.85), hjust = 0.5, color = "gray30"),
      axis.title    = element_text(size = rel(1)),
      axis.text     = element_text(size = rel(0.9)),

      strip.background = element_rect(fill = "gray90", color = "gray60"),
      strip.text       = element_text(size = rel(0.8)),

      legend.title    = element_text(size = rel(1)),
      legend.text     = element_text(size = rel(0.9)),
      legend.position = "bottom",
      legend.box      = "horizontal",

      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),

      plot.margin = margin(10, 10, 10, 10)
    )
}

# ==============================================================================
# DATA LOADING AND PREPROCESSING
# ==============================================================================

log_message("Loading combined ranking outcomes...")

required_cols <- c("original_batch_id", "original_pair_id", "true_population",
                    "true_relationship", "loci_set", "tested_relationship",
                    "tested_population", "combined_LR", "rank", "n_database",
                    "top_200", "top_100", "top_50", "top_10")

outcomes_raw <- fread(input_file)

missing_cols <- setdiff(required_cols, names(outcomes_raw))
if (length(missing_cols) > 0) {
  stop("Input file missing columns: ", paste(missing_cols, collapse = ", "))
}

has_tied_group_size <- "tied_group_size" %in% names(outcomes_raw)
if (!has_tied_group_size) {
  log_message("Note: 'tied_group_size' column not found (older run, before this diagnostic was added) - skipping Figure 3.")
}

log_message(sprintf("Loaded %s rows across %s unique focal pairs",
                     format(nrow(outcomes_raw), big.mark = ","),
                     format(n_distinct(outcomes_raw$original_pair_id), big.mark = ",")))

outcomes <- outcomes_raw %>%
  mutate(
    combined_LR = as.numeric(combined_LR),
    true_relationship = factor(true_relationship,
                               levels = RELATIONSHIP_ORDER,
                               labels = RELATIONSHIP_LABELS),
    tested_relationship = factor(tested_relationship,
                                 levels = RELATIONSHIP_ORDER,
                                 labels = RELATIONSHIP_LABELS),
    loci_set = factor(loci_set,
                      levels = LOCI_ORDER,
                      labels = LOCI_LABELS),
    # Matched-hypothesis indicator (same convention as
    # plots_mismatched_relationship.R's rel_match)
    rel_match = as.character(true_relationship) == as.character(tested_relationship)
  )

# ==============================================================================
# FIGURE 1: PROPORTION OF TRUE RELATIVES RECOVERED IN TOP N
# Gold border highlights bars where the tested hypothesis matches the truth
# ==============================================================================

log_message("Generating Figure 1: Proportion recovered in top N...")

outcomes_long <- outcomes %>%
  pivot_longer(
    cols = c(top_10, top_50, top_100, top_200),
    names_to = "top_n_label",
    values_to = "in_top_n"
  ) %>%
  mutate(
    top_n_label = factor(top_n_label,
                         levels = c("top_10", "top_50", "top_100", "top_200"),
                         labels = c("Top 10", "Top 50", "Top 100", "Top 200"))
  )

prop_summary <- outcomes_long %>%
  group_by(true_relationship, tested_relationship, loci_set, top_n_label) %>%
  summarize(
    prop_in_top_n = mean(in_top_n),
    n_replicates  = n(),
    rel_match     = unique(rel_match),
    .groups = "drop"
  )

fig_prop <- ggplot(
  prop_summary,
  aes(x = top_n_label, y = prop_in_top_n, fill = tested_relationship,
      group = tested_relationship)
) +
  geom_col(
    aes(color = rel_match),
    position = position_dodge(width = 0.8), width = 0.7, linewidth = 0.9
  ) +
  facet_grid(true_relationship ~ loci_set) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_fill_manual(values = relationship_colors, name = "Tested Relationship") +
  scale_color_manual(values = c("TRUE" = "gold", "FALSE" = NA), guide = "none") +
  labs(
    title    = "Focal Ranking Test: Proportion of True Relatives Recovered in Top N",
    subtitle = sprintf("n = %s focal replicates per true relationship | gold border = tested hypothesis matches truth",
                       format(n_distinct(outcomes$original_pair_id), big.mark = ",")),
    x = "Rank Threshold",
    y = "Proportion of True Relatives in Top N"
  ) +
  theme_publication() +
  theme(legend.position = "bottom")

ggsave(
  file.path(output_dir, "ranking_prop_in_top_n.pdf"),
  plot = fig_prop, width = 10, height = 8, units = "in", bg = "white"
)
ggsave(
  file.path(output_dir, "ranking_prop_in_top_n.png"),
  plot = fig_prop, width = 10, height = 8, units = "in", dpi = 300, bg = "white"
)
log_message("Figure 1 saved: ranking_prop_in_top_n.png/.pdf")

# ==============================================================================
# FIGURE 2: DISTRIBUTION OF TRUE RELATIVE RANK (log10)
# Gold shading highlights matched-hypothesis panels
# (same geom_rect technique as plots_mismatched_relationship.R)
# ==============================================================================

log_message("Generating Figure 2: Rank distribution...")

fig_rank <- ggplot(
  outcomes,
  aes(x = tested_relationship, y = rank, fill = tested_relationship)
) +
  geom_rect(
    data = outcomes %>%
      filter(rel_match) %>%
      distinct(true_relationship, tested_relationship, loci_set),
    aes(xmin = as.numeric(tested_relationship) - 0.45,
        xmax = as.numeric(tested_relationship) + 0.45,
        ymin = 1, ymax = max(outcomes$n_database)),
    fill = "gold", alpha = 0.15, inherit.aes = FALSE
  ) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3, alpha = 0.8) +
  facet_grid(true_relationship ~ loci_set) +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = comma_format()
  ) +
  scale_fill_manual(values = relationship_colors, guide = "none") +
  labs(
    title    = "Focal Ranking Test: Distribution of True Relative Rank",
    subtitle = sprintf("n_database = %s candidates per replicate | gold shading = tested hypothesis matches truth",
                       format(unique(outcomes$n_database)[1], big.mark = ",")),
    x = "Tested Relationship",
    y = "Rank of True Relative"
  ) +
  theme_publication() +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(angle = 45, hjust = 1)
  )

ggsave(
  file.path(output_dir, "ranking_rank_distribution.pdf"),
  plot = fig_rank, width = 10, height = 8, units = "in", bg = "white"
)
ggsave(
  file.path(output_dir, "ranking_rank_distribution.png"),
  plot = fig_rank, width = 10, height = 8, units = "in", dpi = 300, bg = "white"
)
log_message("Figure 2 saved: ranking_rank_distribution.png/.pdf")

# ==============================================================================
# FIGURE 3 (OPTIONAL): TIED GROUP SIZE
# Only produced when tied_group_size is present in the input file
# ==============================================================================

if (has_tied_group_size) {
  log_message("Generating Figure 3: Tied group size...")

  fig_tied <- ggplot(
    outcomes,
    aes(x = tested_relationship, y = tied_group_size, fill = tested_relationship)
  ) +
    geom_rect(
      data = outcomes %>%
        filter(rel_match) %>%
        distinct(true_relationship, tested_relationship, loci_set),
      aes(xmin = as.numeric(tested_relationship) - 0.45,
          xmax = as.numeric(tested_relationship) + 0.45,
          ymin = 1, ymax = max(outcomes$n_database)),
      fill = "gold", alpha = 0.15, inherit.aes = FALSE
    ) +
    geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3, alpha = 0.8) +
    facet_grid(true_relationship ~ loci_set) +
    scale_y_log10(
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = comma_format()
    ) +
    scale_fill_manual(values = relationship_colors, guide = "none") +
    labs(
      title    = "Focal Ranking Test: Size of the Tied-Score Group Containing the True Relative",
      subtitle = "Large tied groups (common under a mismatched hypothesis, e.g. combined_LR = 0) mean\nrank is an average over an indistinguishable block, not a singular position",
      x = "Tested Relationship",
      y = "Tied Group Size\n(# candidates sharing the true relative's exact combined LR)"
    ) +
    theme_publication() +
    theme(
      legend.position = "none",
      axis.text.x     = element_text(angle = 45, hjust = 1)
    )

  ggsave(
    file.path(output_dir, "ranking_tied_group_size.pdf"),
    plot = fig_tied, width = 10, height = 8.5, units = "in", bg = "white"
  )
  ggsave(
    file.path(output_dir, "ranking_tied_group_size.png"),
    plot = fig_tied, width = 10, height = 8.5, units = "in", dpi = 300, bg = "white"
  )
  log_message("Figure 3 saved: ranking_tied_group_size.png/.pdf")
}

# ==============================================================================
# SUMMARY STATISTICS TABLE
# ==============================================================================

log_message("Generating summary statistics table...")

summary_stats <- outcomes %>%
  group_by(loci_set, tested_relationship, tested_population, true_relationship, rel_match) %>%
  summarize(
    n_replicates  = n(),
    prop_top_10   = mean(top_10),
    prop_top_50   = mean(top_50),
    prop_top_100  = mean(top_100),
    prop_top_200  = mean(top_200),
    median_rank   = median(rank),
    mean_rank     = mean(rank),
    sd_rank       = sd(rank),
    min_rank      = min(rank),
    max_rank      = max(rank),
    mean_tied_group_size   = if (has_tied_group_size) mean(tied_group_size) else NA_real_,
    median_tied_group_size = if (has_tied_group_size) median(tied_group_size) else NA_real_,
    .groups = "drop"
  ) %>%
  arrange(true_relationship, tested_relationship, loci_set)

summary_filename <- "ranking_summary_statistics.csv"
write_csv(summary_stats, file.path(output_dir, summary_filename))
log_message(sprintf("Summary table saved: %s", summary_filename))

cat("\nMatched-hypothesis recovery (tested_relationship == true_relationship):\n")
print(as.data.frame(
  summary_stats %>%
    filter(rel_match) %>%
    select(true_relationship, loci_set, n_replicates, prop_top_10, prop_top_200, median_rank)
))

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\nOUTPUT FILES:\n")
cat("   - ranking_prop_in_top_n.png/.pdf (Figure 1)\n")
cat("   - ranking_rank_distribution.png/.pdf (Figure 2)\n")
if (has_tied_group_size) {
  cat("   - ranking_tied_group_size.png/.pdf (Figure 3)\n")
}
cat(sprintf("   - %s (summary table)\n", summary_filename))
cat("\nAll files saved to:", output_dir, "\n")
cat("================================================================================\n")

log_message("Analysis complete!")
