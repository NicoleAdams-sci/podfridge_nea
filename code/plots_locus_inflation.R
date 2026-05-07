#!/usr/bin/env Rscript
# =============================================================================
# plots_locus_inflation.R
#
# Purpose: Visualize which loci drive LR inflation under mismatched population
#          allele frequencies, and identify primary driver loci.
#
# Creates:
#   Fig 1:  Ranked locus inflation — AfAm <-> Asian mismatch, full siblings
#           and parent-child. Lollipop chart, loci colored by CODIS panel.
#   Fig 2:  All mismatch directions — heatmap of per-locus median log10 ratio
#           for full siblings across all true x tested population combinations.
#   Fig 3:  Observed heterozygosity for top driver loci by population.
#           Shows biological basis for inflation pattern.
#   Table:  Primary driver loci CSV — loci with median log10 ratio > threshold
#           in at least one key mismatch direction, annotated with panel
#           membership and heterozygosity.
#
# Input:
#   <input_dir>/locus_inflation_summary.csv
#   <input_dir>/locus_heterozygosity_summary.csv
#   data/core_CODIS_loci.csv
#
# Usage:
#   Rscript code/plots_locus_inflation.R <input_dir> [output_dir]
#
#   input_dir   Path to locus_inflation output dir (e.g., output/lr_analysis_20260410/locus_inflation)
#   output_dir  Where to write figures (default: <input_dir>/plots)
#
# Date: 2026-05-07
# =============================================================================

suppressMessages({
  library(tidyverse)
  library(data.table)
  library(scales)
  library(patchwork)
  library(ggrepel)
})

log_message <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time()), msg))

# =============================================================================
# ARGUMENT PARSING
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript code/plots_locus_inflation.R <input_dir> [output_dir]")
}

input_dir  <- args[1]
output_dir <- if (length(args) >= 2) args[2] else file.path(input_dir, "plots")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

log_message(paste("Input directory: ", input_dir))
log_message(paste("Output directory:", output_dir))

# =============================================================================
# CONSTANTS — consistent with other plot scripts
# =============================================================================

RELATIONSHIP_ORDER  <- c("parent_child", "full_siblings", "half_siblings",
                         "cousins", "second_cousins", "unrelated")
RELATIONSHIP_LABELS <- c("Parent-Child", "Full Siblings", "Half Siblings",
                         "Cousins", "Second Cousins", "Unrelated")
names(RELATIONSHIP_LABELS) <- RELATIONSHIP_ORDER

# Okabe-Ito population colors (matches plots_matched_publication.R)
POPULATION_COLORS <- c(
  "AfAm"     = "#0072B2",
  "Asian"    = "#009E73",
  "Cauc"     = "#56B4E9",
  "Hispanic" = "#CC79A7"
)
POPULATION_LABELS <- c(
  "AfAm"     = "African American",
  "Asian"    = "Asian",
  "Cauc"     = "Caucasian",
  "Hispanic" = "Hispanic"
)

# CODIS panel colors — from least to most extended
PANEL_ORDER <- c("Core 13", "Identifiler 15", "Expanded 20",
                 "Supplementary", "Autosomal 29 only")
PANEL_COLORS <- c(
  "Core 13"          = "#E69F00",
  "Identifiler 15"   = "#56B4E9",
  "Expanded 20"      = "#009E73",
  "Supplementary"    = "#CC79A7",
  "Autosomal 29 only"= "#D55E00"
)

# Inflation threshold for "primary driver" classification
DRIVER_THRESHOLD <- 0.05  # median log10 ratio

# Publication theme — matches other scripts
theme_publication <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title       = element_text(size = rel(1.2), hjust = 0.5),
      plot.subtitle    = element_text(size = rel(1.0), hjust = 0.5, color = "gray30"),
      axis.title       = element_text(size = rel(1.0)),
      axis.text        = element_text(size = rel(0.9)),
      strip.background = element_rect(fill = "gray90", color = "gray60"),
      strip.text       = element_text(size = rel(0.85)),
      legend.title     = element_text(size = rel(1.0)),
      legend.text      = element_text(size = rel(0.9)),
      legend.position  = "bottom",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      plot.margin      = margin(10, 10, 10, 10)
    )
}

# =============================================================================
# DATA LOADING
# =============================================================================

log_message("Loading data...")

inflation_file <- file.path(input_dir, "locus_inflation_summary.csv")
heteroz_file   <- file.path(input_dir, "locus_heterozygosity_summary.csv")
codis_file     <- "data/core_CODIS_loci.csv"

for (f in c(inflation_file, heteroz_file, codis_file)) {
  if (!file.exists(f)) stop(paste("File not found:", f))
}

inflation <- fread(inflation_file)
heteroz   <- fread(heteroz_file)
codis     <- fread(codis_file)

log_message(sprintf("Loaded %d inflation rows, %d loci in heterozygosity table.",
                    nrow(inflation), nrow(heteroz)))

# =============================================================================
# ASSIGN CODIS PANEL MEMBERSHIP
# "First appearance" panel — most informative for forensic context
# Loci not in core_CODIS_loci.csv are autosomal_29 only
# =============================================================================

panel_lookup <- codis %>%
  mutate(
    panel = case_when(
      core_13        == 1 ~ "Core 13",
      identifiler_15 == 1 ~ "Identifiler 15",
      expanded_20    == 1 ~ "Expanded 20",
      supplementary  == 1 ~ "Supplementary",
      TRUE                ~ "Autosomal 29 only"
    ),
    panel = factor(panel, levels = PANEL_ORDER)
  ) %>%
  select(locus, panel)

# Loci only in autosomal_29 (not in codis CSV at all)
all_loci     <- unique(inflation$locus)
missing_loci <- setdiff(all_loci, codis$locus)
if (length(missing_loci) > 0) {
  panel_lookup <- bind_rows(
    panel_lookup,
    data.frame(locus = missing_loci,
               panel = factor("Autosomal 29 only", levels = PANEL_ORDER))
  )
}

inflation <- inflation %>% left_join(panel_lookup, by = "locus")
heteroz   <- heteroz   %>% left_join(panel_lookup, by = "locus")

# Apply factor levels
inflation <- inflation %>%
  mutate(
    known_relationship = factor(known_relationship,
                                levels = RELATIONSHIP_ORDER,
                                labels = RELATIONSHIP_LABELS),
    panel = factor(panel, levels = PANEL_ORDER)
  )

# =============================================================================
# FIGURE 1: RANKED LOCUS INFLATION — AfAm <-> Asian
# Lollipop chart for full siblings and parent-child
# =============================================================================

log_message("Generating Figure 1: Ranked locus inflation (AfAm <-> Asian)...")

fig1_data <- inflation %>%
  filter(
    ((population == "AfAm" & tested_population == "Asian") |
     (population == "Asian" & tested_population == "AfAm")),
    known_relationship %in% c("Full Siblings", "Parent-Child")
  ) %>%
  mutate(
    direction = if_else(population == "AfAm",
                        "True AfAm, tested with\nAsian frequencies",
                        "True Asian, tested with\nAfAm frequencies")
  )

# Rank loci within each direction x relationship panel by median inflation
# Use full_siblings AfAm->Asian as the canonical rank order
canonical_rank <- fig1_data %>%
  filter(direction == "True AfAm, tested with\nAsian frequencies",
         known_relationship == "Full Siblings") %>%
  arrange(desc(median_log10_ratio)) %>%
  mutate(locus_ordered = factor(locus, levels = rev(locus)))

locus_order <- levels(canonical_rank$locus_ordered)

fig1_data <- fig1_data %>%
  mutate(locus_ordered = factor(locus, levels = locus_order))

fig1 <- ggplot(fig1_data,
               aes(x = median_log10_ratio, y = locus_ordered, color = panel)) +

  # Zero reference line
  geom_vline(xintercept = 0, color = "gray60", linewidth = 0.4, linetype = "dashed") +

  # Driver threshold line
  geom_vline(xintercept = DRIVER_THRESHOLD, color = "gray40", linewidth = 0.4,
             linetype = "dotted") +

  # CI segment
  geom_segment(aes(x = ci_lower, xend = ci_upper,
                   y = locus_ordered, yend = locus_ordered),
               alpha = 0.3, linewidth = 0.6) +

  # Point
  geom_point(size = 2.5) +

  facet_grid(direction ~ known_relationship) +

  scale_color_manual(values = PANEL_COLORS, name = "CODIS Panel",
                     drop = FALSE) +

  labs(
    title    = "Per-Locus LR Inflation Under African American \u2194 Asian Frequency Mismatch",
    subtitle = paste0("Median log\u2081\u2080(LR_wrong / LR_correct) per locus | True positives only\n",
                      "Dotted line = driver threshold (", DRIVER_THRESHOLD, ")"),
    x        = "Median log\u2081\u2080(LR_wrong / LR_correct)",
    y        = NULL
  ) +

  theme_publication() +
  theme(
    legend.position = "right",
    axis.text.y     = element_text(size = rel(0.75))
  )

fig1_file <- "locus_inflation_ranked.png"
ggsave(file.path(output_dir, fig1_file),
       plot = fig1, width = 12, height = 10, dpi = 300, bg = "white")
log_message(paste("Figure 1 saved:", fig1_file))

# =============================================================================
# FIGURE 2: HEATMAP — all mismatch directions, full siblings
# =============================================================================

log_message("Generating Figure 2: All-direction mismatch heatmap...")

fig2_data <- inflation %>%
  filter(
    known_relationship == "Full Siblings",
    population != tested_population   # mismatched only
  ) %>%
  mutate(
    true_pop   = factor(population,        levels = names(POPULATION_LABELS),
                        labels = POPULATION_LABELS),
    tested_pop = factor(tested_population, levels = names(POPULATION_LABELS),
                        labels = POPULATION_LABELS),
    # Use same factor levels as Fig 1 so F13B appears at top in both figures
    locus_ordered = factor(locus, levels = locus_order)
  )

fig2 <- ggplot(fig2_data,
               aes(x = tested_pop, y = locus_ordered,
                   fill = median_log10_ratio)) +

  geom_tile(color = "white", linewidth = 0.4) +

  facet_wrap(~ true_pop, nrow = 1) +

  scale_fill_gradientn(
    colors = c("white", "#c6dbef", "#6baed6", "#2166ac", "#08306b"),
    limits = c(0, NA),
    oob    = scales::squish,
    name   = "Median\nlog\u2081\u2080(LR_wrong /\nLR_correct)"
  ) +

  labs(
    title    = "Per-Locus LR Inflation Across All Population Mismatch Directions",
    subtitle = "Full siblings | Each cell = one locus x tested frequency database combination",
    x        = "Tested Population Frequencies",
    y        = NULL
  ) +

  theme_publication() +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = rel(0.8)),
    axis.text.y     = element_text(size  = rel(0.7)),
    panel.grid      = element_blank(),
    legend.position = "right"
  )

fig2_file <- "locus_inflation_heatmap.png"
ggsave(file.path(output_dir, fig2_file),
       plot = fig2, width = 14, height = 10, dpi = 300, bg = "white")
log_message(paste("Figure 2 saved:", fig2_file))

# =============================================================================
# FIGURE 3: HETEROZYGOSITY DIFFERENCE vs. LR INFLATION — ALL 29 LOCI
#
# Scatter plot: x = observed heterozygosity difference (AfAm - Asian),
#               y = median log10 inflation (AfAm->Asian, full siblings)
# All 29 loci included. Driver loci labeled. Trend line added.
# Honest framing: shows where heterozygosity explains inflation and where
# it does not — loci with negative x (Asian more heterozygous) but positive
# y (still inflating) are visible and not hidden.
# =============================================================================

log_message("Generating Figure 3: Heterozygosity difference vs. LR inflation (all loci)...")

# Identify driver loci for labeling (used in table too)
driver_loci <- inflation %>%
  filter(
    known_relationship == "Full Siblings",
    ((population == "AfAm" & tested_population == "Asian") |
     (population == "Asian" & tested_population == "AfAm"))
  ) %>%
  group_by(locus) %>%
  summarise(max_median = max(median_log10_ratio), .groups = "drop") %>%
  filter(max_median >= DRIVER_THRESHOLD) %>%
  arrange(desc(max_median)) %>%
  pull(locus)

log_message(sprintf("Identified %d primary driver loci at threshold %.2f.",
                    length(driver_loci), DRIVER_THRESHOLD))

# Build scatter data: all 29 loci
# x = AfAm - Asian heterozygosity difference (from heteroz table)
# y = median log10 inflation AfAm->Asian, full siblings
fig3_data <- inflation %>%
  filter(population == "AfAm", tested_population == "Asian",
         known_relationship == "Full Siblings") %>%
  select(locus, panel, median_log10_ratio) %>%
  left_join(
    heteroz %>% select(locus, het_AfAm = AfAm, het_Asian = Asian,
                       het_diff = AfAm_minus_Asian),
    by = "locus"
  ) %>%
  mutate(is_driver = locus %in% driver_loci)

fig3 <- ggplot(fig3_data,
               aes(x = het_diff, y = median_log10_ratio, color = panel)) +

  # Reference lines
  geom_hline(yintercept = 0,              color = "gray70", linewidth = 0.4) +
  geom_hline(yintercept = DRIVER_THRESHOLD, color = "gray40", linewidth = 0.4,
             linetype = "dotted") +
  geom_vline(xintercept = 0,              color = "gray70", linewidth = 0.4) +

  # Trend line across all loci — shows overall relationship
  geom_smooth(aes(group = 1), method = "lm", se = TRUE,
              color = "gray50", fill = "gray85", linewidth = 0.7,
              show.legend = FALSE) +

  # All loci
  geom_point(size = 2.5, alpha = 0.85) +

  # Label driver loci only
  ggrepel::geom_text_repel(
    data    = filter(fig3_data, is_driver),
    aes(label = locus),
    size    = 3.2,
    color   = "gray20",
    box.padding     = 0.4,
    point.padding   = 0.3,
    min.segment.len = 0.2,
    show.legend     = FALSE
  ) +

  scale_color_manual(values = PANEL_COLORS, name = "CODIS Panel", drop = FALSE) +

  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +

  labs(
    title    = "Heterozygosity Difference vs. LR Inflation: All 29 Loci",
    subtitle = paste0(
      "x-axis: observed heterozygosity (AfAm \u2212 Asian) | ",
      "y-axis: median log\u2081\u2080(LR_wrong / LR_correct), AfAm\u2192Asian, full siblings\n",
      "Labeled loci exceed driver threshold (",  DRIVER_THRESHOLD, "). ",
      "Loci with negative x show Asian is more heterozygous yet some still inflate."
    ),
    x = "Heterozygosity Difference (AfAm \u2212 Asian)",
    y = "Median log\u2081\u2080(LR_wrong / LR_correct)"
  ) +

  theme_publication() +
  theme(legend.position = "right")

fig3_file <- "locus_inflation_heterozygosity.png"
ggsave(file.path(output_dir, fig3_file),
       plot = fig3, width = 10, height = 7, dpi = 300, bg = "white")
log_message(paste("Figure 3 saved:", fig3_file))

# =============================================================================
# TABLE: PRIMARY DRIVER LOCI
# One row per locus, annotated with panel, inflation (both directions),
# heterozygosity (AfAm, Asian, difference)
# =============================================================================

log_message("Building primary driver loci table...")

# Pull key inflation columns: AfAm->Asian and Asian->AfAm, full siblings
inflation_afam_asian <- inflation %>%
  filter(population == "AfAm", tested_population == "Asian",
         known_relationship == "Full Siblings") %>%
  select(locus, panel,
         fs_AfAm_Asian_median  = median_log10_ratio,
         fs_AfAm_Asian_ci_lower = ci_lower,
         fs_AfAm_Asian_ci_upper = ci_upper)

inflation_asian_afam <- inflation %>%
  filter(population == "Asian", tested_population == "AfAm",
         known_relationship == "Full Siblings") %>%
  select(locus,
         fs_Asian_AfAm_median  = median_log10_ratio,
         fs_Asian_AfAm_ci_lower = ci_lower,
         fs_Asian_AfAm_ci_upper = ci_upper)

# Parent-child AfAm->Asian for completeness
inflation_pc <- inflation %>%
  filter(population == "AfAm", tested_population == "Asian",
         known_relationship == "Parent-Child") %>%
  select(locus,
         pc_AfAm_Asian_median = median_log10_ratio)

# Heterozygosity columns
heteroz_cols <- heteroz %>%
  select(locus, het_AfAm = AfAm, het_Asian = Asian,
         het_Cauc = Cauc, het_Hispanic = Hispanic,
         het_AfAm_minus_Asian = AfAm_minus_Asian)

# Combine
driver_table <- inflation_afam_asian %>%
  left_join(inflation_asian_afam, by = "locus") %>%
  left_join(inflation_pc,         by = "locus") %>%
  left_join(heteroz_cols,         by = "locus") %>%
  mutate(
    is_primary_driver = (fs_AfAm_Asian_median >= DRIVER_THRESHOLD |
                         fs_Asian_AfAm_median >= DRIVER_THRESHOLD),
    driver_threshold  = DRIVER_THRESHOLD
  ) %>%
  arrange(desc(fs_AfAm_Asian_median))

# Round numeric columns for readability
driver_table <- driver_table %>%
  mutate(across(where(is.numeric) & !matches("threshold"), ~ round(.x, 4)))

table_file <- file.path(output_dir, "primary_driver_loci_table.csv")
fwrite(driver_table, table_file)
log_message(sprintf("Driver loci table saved: %s (%d loci, %d flagged as primary drivers)",
                    table_file,
                    nrow(driver_table),
                    sum(driver_table$is_primary_driver)))

# =============================================================================
# CONSOLE SUMMARY
# =============================================================================

cat("\n")
cat("=============================================================================\n")
cat("PRIMARY DRIVER LOCI (median log10 ratio >= ", DRIVER_THRESHOLD,
    " in AfAm<->Asian full siblings)\n", sep = "")
cat("=============================================================================\n")

primary <- driver_table %>%
  filter(is_primary_driver) %>%
  select(locus, panel,
         fs_AfAm_Asian_median, fs_Asian_AfAm_median,
         het_AfAm, het_Asian, het_AfAm_minus_Asian)

print(as.data.frame(primary))

cat("\n")
cat("OUTPUT FILES:\n")
for (f in c(fig1_file, fig2_file, fig3_file, "primary_driver_loci_table.csv")) {
  cat(sprintf("  %s\n", file.path(output_dir, f)))
}
cat("=============================================================================\n")

log_message("Done.")
