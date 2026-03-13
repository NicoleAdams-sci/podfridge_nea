#!/usr/bin/env Rscript
# ==============================================================================
# Publication-Ready Plots for Population & Relationship Matched LR Analysis
# ==============================================================================
# Creates:
#   1. Main text figure: Violin plots by relationship type (all populations)
#   2. Supplement figure: Box plots by population (showing consistency)
#
# Date: 2026-03-06
# ==============================================================================

# Load Required Libraries ----
suppressMessages(suppressWarnings({
  library(tidyverse)
  library(data.table)
  library(scales)
  library(patchwork)  # For combining plots if needed
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
COMBINED_LR_MATCH_FILE <- "combined_LR_match.csv"

# Factor levels for consistent ordering
relationship_order <- c("parent_child", "full_siblings", "half_siblings", 
                        "cousins", "second_cousins", "unrelated")
relationship_labels <- c("Parent-Child", "Full Siblings", "Half Siblings", 
                         "Cousins", "Second Cousins", "Unrelated")
names(relationship_labels) <- relationship_order

loci_set_order <- c("core_13", "identifiler_15", "expanded_20", 
                    "supplementary", "autosomal_29")
loci_set_labels <- c("Core 13", "Identifiler 15", "Expanded 20", 
                     "Supplementary", "Autosomal 29")
names(loci_set_labels) <- loci_set_order

# Color palettes ----
# For relationship types - colorblind-friendly palette
relationship_colors <- c(
  "Parent-Child" = "#D55E00",      # Vermillion
  "Full Siblings" = "#E69F00",     # Orange  
  "Half Siblings" = "#56B4E9",     # Sky blue
  "Cousins" = "#009E73",           # Bluish green
  "Second Cousins" = "#CC79A7",    # Reddish purple
  "Unrelated" = "#999999"          # Gray
)

# For populations
population_colors <- c(
  "AfAm" = "#E41A1C", 
  "Asian" = "#377EB8", 
  "Cauc" = "#4DAF4A", 
  "Hispanic" = "#984EA3", 
  "all" = "#FF7F00"
)

# Data Loading ----
log_message("Loading matched LR data...")

raw_lrs_path <- file.path(input_dir, COMBINED_LR_MATCH_FILE)
combined_lrs_match <- fread(raw_lrs_path) %>%
  as_tibble() %>%
  # Filter for population and relationship matched only
  filter(is_correct_pop == TRUE, known_relationship == tested_relationship) %>%
  # Apply factor labels
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
    population = factor(population, levels = c("AfAm", "Asian", "Cauc", "Hispanic", "all")),
    log_LR = log10(combined_LR)
  )

log_message(paste("Loaded", nrow(combined_lrs_match), "matched LR observations"))

# ==============================================================================
# MAIN TEXT FIGURE: All Populations Combined
# ==============================================================================

log_message("Creating main text figure...")

# Filter for combined population data
main_data <- combined_lrs_match %>%
  filter(population == "all")

# Create violin plot with clean publication styling
fig_main <- ggplot(main_data, aes(x = relationship, y = log_LR, fill = relationship)) +
  # Violin plots showing full distribution
  geom_violin(
    alpha = 0.7,
    draw_quantiles = c(0.25, 0.5, 0.75),
    scale = "width",
    trim = TRUE
  ) +
  # Add mean as a point
  stat_summary(
    fun = mean,
    geom = "point",
    size = 2,
    color = "black",
    shape = 18  # Diamond shape
  ) +
  # Facet by loci set
  facet_wrap(~ loci_set, ncol = 5, scales = "fixed") +
  # Apply color scheme
  scale_fill_manual(values = relationship_colors) +
  # Format y-axis as 10^x
  scale_y_continuous(
    breaks = seq(0, 40, 10),
    labels = trans_format("identity", math_format(10^.x)),
    limits = c(-5, 40)
  ) +
  # Labels
  labs(
    x = NULL,
    y = expression(paste("Combined Likelihood Ratio (", log[10], " scale)")),
    fill = "Relationship Type"
  ) +
  # Clean theme
  theme_bw(base_size = 11) +
  theme(
    # Remove x-axis text since legend shows relationships
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    # Format axis text
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 11, face = "bold"),
    # Format facet labels
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "gray90", color = "gray60"),
    # Legend at bottom
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    # Panel formatting
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.8, "lines"),
    # Plot margins
    plot.margin = margin(10, 10, 10, 10)
  ) +
  # Format legend
  guides(fill = guide_legend(nrow = 1, override.aes = list(alpha = 0.8)))

# Save main figure
ggsave(
  file.path(output_dir, "matched_main_lr_distributions.pdf"),
  plot = fig_main,
  width = 12,
  height = 5,
  units = "in",
  bg = "white"
)

ggsave(
  file.path(output_dir, "matched_main_lr_distributions.png"),
  plot = fig_main,
  width = 12,
  height = 5,
  units = "in",
  dpi = 300,
  bg = "white"
)

log_message("Main text figure saved")

# ==============================================================================
# SUPPLEMENT FIGURE: By Population
# ==============================================================================

log_message("Creating supplement figure...")

# Use all populations including "all" (which uses combined allele frequencies)
supp_data <- combined_lrs_match

# Create box plot showing consistency across populations
fig_supp <- ggplot(supp_data, aes(x = relationship, y = log_LR, fill = population)) +
  # Box plots with smaller outlier points
  geom_boxplot(
    position = position_dodge(width = 0.85),
    alpha = 0.75,
    outlier.size = 0.3,
    outlier.alpha = 0.3,
    linewidth = 0.4
  ) +
  # Facet by loci set
  facet_wrap(~ loci_set, ncol = 3, scales = "fixed") +
  # Apply population colors
  scale_fill_manual(
    values = population_colors,
    labels = c("African American", "Asian", "Caucasian", "Hispanic", "All Populations")
  ) +
  # Format y-axis
  scale_y_continuous(
    breaks = seq(0, 40, 10),
    labels = trans_format("identity", math_format(10^.x)),
    limits = c(-5, 40)
  ) +
  # Labels
  labs(
    x = "Relationship Type",
    y = expression(paste("Combined Likelihood Ratio (", log[10], " scale)")),
    fill = "Population"
  ) +
  # Clean theme
  theme_bw(base_size = 11) +
  theme(
    # Format axis text
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 11, face = "bold"),
    # Format facet labels
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "gray90", color = "gray60"),
    # Legend at bottom
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    # Panel formatting
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1, "lines"),
    # Plot margins
    plot.margin = margin(10, 10, 10, 10)
  ) +
  # Format legend
  guides(fill = guide_legend(nrow = 1, override.aes = list(alpha = 0.8)))

# Save supplement figure
ggsave(
  file.path(output_dir, "matched_supp_lr_by_population.pdf"),
  plot = fig_supp,
  width = 12,
  height = 10,
  units = "in",
  bg = "white"
)

ggsave(
  file.path(output_dir, "matched_supp_lr_by_population.png"),
  plot = fig_supp,
  width = 12,
  height = 10,
  units = "in",
  dpi = 300,
  bg = "white"
)

log_message("Supplement figure saved")

# ==============================================================================
# Summary Statistics for Results Section
# ==============================================================================

log_message("Calculating summary statistics for manuscript...")

# Calculate key statistics for the "all" population
summary_stats <- main_data %>%
  group_by(relationship, loci_set) %>%
  summarise(
    mean_logLR = mean(log_LR),
    median_logLR = median(log_LR),
    sd_logLR = sd(log_LR),
    q25 = quantile(log_LR, 0.25),
    q75 = quantile(log_LR, 0.75),
    .groups = 'drop'
  ) %>%
  arrange(loci_set, desc(mean_logLR))

# Save summary statistics
write_csv(
  summary_stats,
  file.path(output_dir, "matched_summary_statistics_for_ms.csv")
)

# Print key findings to console
log_message("\n=== Key Findings for Results Section ===")

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

cat("\nUnrelated mean log10(LR) by loci set (should be ~0):\n")
print(unrel_stats, n = Inf)

# Calculate separation between adjacent relationship types
# (e.g., Full Sibs vs Half Sibs)
separation_stats <- summary_stats %>%
  arrange(loci_set, relationship) %>%
  group_by(loci_set) %>%
  mutate(
    separation_from_next = mean_logLR - lead(mean_logLR)
  ) %>%
  filter(!is.na(separation_from_next)) %>%
  select(loci_set, relationship, mean_logLR, separation_from_next)

cat("\nSeparation between adjacent relationship types:\n")
print(separation_stats, n = Inf)

# Population consistency check
pop_consistency <- combined_lrs_match %>%
  filter(population != "all") %>%
  group_by(relationship, loci_set, population) %>%
  summarise(mean_logLR = mean(log_LR), .groups = 'drop') %>%
  group_by(relationship, loci_set) %>%
  summarise(
    mean_across_pops = mean(mean_logLR),
    sd_across_pops = sd(mean_logLR),
    cv_across_pops = sd(mean_logLR) / mean(mean_logLR),
    .groups = 'drop'
  ) %>%
  arrange(loci_set, desc(mean_across_pops))

cat("\nPopulation consistency (low CV = consistent across populations):\n")
print(pop_consistency, n = Inf)

write_csv(
  pop_consistency,
  file.path(output_dir, "population_consistency_check.csv")
)

log_message("\n=== Script Complete ===")
log_message(paste("Figures saved to:", output_dir))
