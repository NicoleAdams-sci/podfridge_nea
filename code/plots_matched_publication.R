#!/usr/bin/env Rscript
# ==============================================================================
# Publication-Ready Plots for Population & Relationship Matched LR Analysis
# ==============================================================================
# Creates:
#   1. Main text figure: Violin plots by relationship type (all populations)
#   2. Supplement figure: Box plots by population (showing consistency)
#   3. Descriptive summary statistics CSV for results section
#
# Statistical tests have been moved to run_statistical_tests.R
#
# Date: 2026-03-06
# ==============================================================================

# Load Required Libraries ----
suppressMessages(suppressWarnings({
  library(tidyverse)
  library(data.table)
  library(scales)
  library(patchwork)
}))

# Helper Functions ----
log_message <- function(message) {
  cat(paste0("[", Sys.time(), "] ", message, "\n"))
}

# Argument Parsing ----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript plots_matched_publication.R <input_dir> [output_dir]")
}

input_dir <- args[1]
log_message(paste("Input directory:", input_dir))

if (length(args) >= 2) {
  output_dir <- args[2]
} else {
  output_dir <- file.path(input_dir, "publication_plots")
}
log_message(paste("Output directory:", output_dir))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Define Constants ----
COMBINED_LR_MATCH_FILE <- "combined_LR_match.csv.gz"

# Factor levels for consistent ordering
relationship_order  <- c("parent_child", "full_siblings", "half_siblings",
                         "cousins", "second_cousins", "unrelated")
relationship_labels <- c("Parent-Child", "Full Siblings", "Half Siblings",
                         "Cousins", "Second Cousins", "Unrelated")
names(relationship_labels) <- relationship_order

loci_set_order  <- c("core_13", "identifiler_15", "expanded_20",
                     "supplementary", "autosomal_29")
loci_set_labels <- c("Core 13", "Identifiler 15", "Expanded 20",
                     "Supplementary", "Autosomal 29")
names(loci_set_labels) <- loci_set_order

# Approximate locus counts per panel (used for Spearman correlation)
# Supplementary is approximate; use rank order if exact count is uncertain
loci_counts <- c(
  "Core 13"       = 13,
  "Identifiler 15" = 15,
  "Expanded 20"   = 20,
  "Supplementary" = 23,
  "Autosomal 29"  = 29
)

# ==============================================================================
# COLOR PALETTES — Shared across scripts
# ==============================================================================
# Both palettes draw from the Okabe-Ito colorblind-safe set.
# Relationship and population colors are intentionally non-overlapping so
# figures using both can be placed side by side without confusion.
#
# Okabe-Ito full palette (8 colors):
#   #E69F00  #56B4E9  #009E73  #F0E442  #0072B2  #D55E00  #CC79A7  #000000
#
# Relationship colors use: #D55E00 #E69F00 #56B4E9 #009E73 #CC79A7 #999999
# Population colors use:   #0072B2 #009E73 #56B4E9 #CC79A7 #999999
# (overlap in population/relationship is acceptable — they never share a plot)

# Relationship type colors (Okabe-Ito — warm to cool, most to least related)
relationship_colors <- c(
  "Parent-Child"   = "#D55E00",   # Vermillion
  "Full Siblings"  = "#E69F00",   # Orange
  "Half Siblings"  = "#56B4E9",   # Sky blue
  "Cousins"        = "#009E73",   # Bluish green
  "Second Cousins" = "#CC79A7",   # Reddish purple
  "Unrelated"      = "#999999"    # Gray
)

# Population colors (Okabe-Ito subset — distinct from relationship palette)
population_colors <- c(
  "AfAm"     = "#0072B2",   # Deep blue
  "Asian"    = "#009E73",   # Bluish green
  "Cauc"     = "#56B4E9",   # Sky blue
  "Hispanic" = "#CC79A7",   # Reddish purple
  "all"      = "#999999"    # Gray (combined)
)

population_labels <- c(
  "AfAm"     = "African American",
  "Asian"    = "Asian",
  "Cauc"     = "Caucasian",
  "Hispanic" = "Hispanic",
  "all"      = "All Populations"
)

# ==============================================================================
# Data Loading ----
# ==============================================================================
log_message("Loading matched LR data...")

raw_lrs_path <- file.path(input_dir, COMBINED_LR_MATCH_FILE)
combined_lrs_match <- read_csv(raw_lrs_path, show_col_types = FALSE) %>%
  filter(is_correct_pop == TRUE, known_relationship == tested_relationship) %>%
  mutate(
    combined_LR = as.numeric(combined_LR),
    relationship = factor(
      known_relationship,
      levels = relationship_order,
      labels = relationship_labels
    ),
    loci_set = factor(
      loci_set,
      levels = loci_set_order,
      labels = loci_set_labels
    ),
    population = factor(
      population,
      levels = c("AfAm", "Asian", "Cauc", "Hispanic", "all"),
      labels = c("AfAm", "Asian", "Cauc", "Hispanic", "all")
    ),
    log_LR = log10(combined_LR)
  )

log_message(paste("Loaded", nrow(combined_lrs_match), "matched LR observations"))

# ==============================================================================
# MAIN TEXT FIGURE: All Populations Combined
# ==============================================================================
log_message("Creating main text figure...")

main_data <- combined_lrs_match %>%
  filter(population == "all")

fig_main <- ggplot(main_data, aes(x = relationship, y = log_LR, fill = relationship)) +
  geom_violin(
    alpha = 0.7,
    draw_quantiles = c(0.25, 0.5, 0.75),
    scale = "width",
    trim = TRUE
  ) +
  stat_summary(
    fun  = mean,
    geom = "point",
    size = 2.5,
    color = "black",
    shape = 18
  ) +
  facet_wrap(~ loci_set, ncol = 5, scales = "fixed") +
  scale_fill_manual(values = relationship_colors) +
  scale_y_continuous(
    breaks = seq(0, 40, 10),
    labels = trans_format("identity", math_format(10^.x)),
    limits = c(-5, 40)
  ) +
  labs(
    x    = NULL,
    y    = expression(paste("Combined Likelihood Ratio (", log[10], " scale)")),
    fill = "Relationship Type"
  ) +
  theme_bw(base_size = 14) +   # increased from 11
  theme(
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    axis.text.y        = element_text(size = 13),
    axis.title.y       = element_text(size = 14),
    strip.text         = element_text(size = 13),
    strip.background   = element_rect(fill = "gray90", color = "gray60"),
    legend.position    = "bottom",
    legend.title       = element_text(size = 13),
    legend.text        = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.spacing      = unit(0.8, "lines"),
    plot.margin        = margin(10, 10, 10, 10)
  ) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(alpha = 0.8)))

ggsave(
  file.path(output_dir, "matched_main_lr_distributions.pdf"),
  plot = fig_main, width = 12, height = 5, units = "in", bg = "white"
)
ggsave(
  file.path(output_dir, "matched_main_lr_distributions.png"),
  plot = fig_main, width = 12, height = 5, units = "in", dpi = 300, bg = "white"
)
log_message("Main text figure saved")

# ==============================================================================
# SUPPLEMENT FIGURE: By Population
# ==============================================================================
log_message("Creating supplement figure...")

supp_data <- combined_lrs_match

fig_supp <- ggplot(supp_data, aes(x = relationship, y = log_LR, fill = population)) +
  geom_boxplot(
    position     = position_dodge(width = 0.85),
    alpha        = 0.75,
    outlier.size  = 0.3,
    outlier.alpha = 0.3,
    linewidth    = 0.4
  ) +
  facet_wrap(~ loci_set, ncol = 3, scales = "fixed") +
  scale_fill_manual(
    values = population_colors,
    labels = population_labels
  ) +
  scale_y_continuous(
    breaks = seq(0, 40, 10),
    labels = trans_format("identity", math_format(10^.x)),
    limits = c(-5, 40)
  ) +
  labs(
    x    = "Relationship Type",
    y    = expression(paste("Combined Likelihood Ratio (", log[10], " scale)")),
    fill = "Population"
  ) +
  theme_bw(base_size = 14) +   # increased from 11
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y      = element_text(size = 13),
    axis.title       = element_text(size = 14),
    strip.text       = element_text(size = 13),
    strip.background = element_rect(fill = "gray90", color = "gray60"),
    legend.position  = "bottom",
    legend.title     = element_text(size = 13),
    legend.text      = element_text(size = 12),
    panel.grid.minor = element_blank(),
    panel.spacing    = unit(1, "lines"),
    plot.margin      = margin(10, 10, 10, 10)
  ) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(alpha = 0.8)))

ggsave(
  file.path(output_dir, "matched_supp_lr_by_population.pdf"),
  plot = fig_supp, width = 12, height = 10, units = "in", bg = "white"
)
ggsave(
  file.path(output_dir, "matched_supp_lr_by_population.png"),
  plot = fig_supp, width = 12, height = 10, units = "in", dpi = 300, bg = "white"
)
log_message("Supplement figure saved")

# Descriptive Summary Statistics
# ==============================================================================
log_message("Calculating descriptive summary statistics...")

summary_stats <- main_data %>%
  group_by(relationship, loci_set) %>%
  summarise(
    mean_logLR   = mean(log_LR),
    median_logLR = median(log_LR),
    sd_logLR     = sd(log_LR),
    q25          = quantile(log_LR, 0.25),
    q75          = quantile(log_LR, 0.75),
    .groups = "drop"
  ) %>%
  arrange(loci_set, desc(mean_logLR))

write_csv(summary_stats,
          file.path(output_dir, "matched_summary_statistics_for_ms.csv"))

# Parent-child LRs across loci sets
pc_stats <- summary_stats %>%
  filter(relationship == "Parent-Child") %>%
  select(loci_set, mean_logLR, median_logLR)

cat("\nParent-Child mean log10(LR) by loci set:\n")
print(pc_stats, n = Inf)

# Unrelated LRs (should be near 0)
unrel_stats <- summary_stats %>%
  filter(relationship == "Unrelated") %>%
  select(loci_set, mean_logLR, sd_logLR)

cat("\nUnrelated mean log10(LR) by loci set (expected ~0):\n")
print(unrel_stats, n = Inf)


log_message("=== Script Complete ===")
log_message(paste("All outputs saved to:", output_dir))
