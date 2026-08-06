#!/usr/bin/env Rscript

# =============================================================================
# Manuscript Figures: LR Classification Performance
# =============================================================================
# Produces three figures:
#
#   FIGURE 1 (Main Text):
#     Classification summary using the 0.01% FPR threshold only, contrasting
#     two loci panels as stacked A/B sub-panels: Core 13 (A) and Expanded 20 (B).
#     Clean, single-threshold view for the main text argument, showing how the
#     larger panel improves discrimination.
#     Output: cutoff_classification_0.01fpr_core13_expanded20.png
#
#   FIGURE 2 (Supplement):
#     Classification summary at 29 autosomal loci across all three FPR thresholds
#     (1%, 0.1%, 0.01%), showing how performance changes with threshold stringency.
#     Output: cutoff_classification_fpr_29loci.png
#
#   FIGURE 3 (Supplement):
#     Heatmap of false positive rates across all loci panels and all four LR
#     thresholds, stratified by population and tested hypothesis.
#     Output: cutoff_supp_heatmap_fp_rates.png
#
# Input:  <input_dir>/proportions_with_classification.csv
#         (produced by prepare_combined_lr_intermediates.R)
#
# Usage:
#   Rscript code/plots_cutoffs_publication.R <input_dir> [output_dir]
#
#   input_dir   Full path to data directory (e.g., output/lr_analysis_20260130)
#   output_dir  Where to write figures (default: <input_dir>/plots_cutoffs)
# =============================================================================

suppressMessages(suppressWarnings({
  library(tidyverse)
  library(data.table)
  library(scales)
  library(patchwork)
}))

log_message <- function(msg) cat(paste0("[", Sys.time(), "] ", msg, "\n"))

# --- Argument parsing ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript code/plots_cutoffs_publication.R <input_dir> [output_dir]")

input_dir  <- args[1]
output_dir <- if (length(args) >= 2) args[2] else file.path(input_dir, "plots_cutoffs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

log_message(paste("Input directory: ", input_dir))
log_message(paste("Output directory:", output_dir))


# =============================================================================
# SECTION 1: SHARED LABELS, FACTORS, AND THEME
# =============================================================================

# Factor orders and human-readable labels used consistently across both figures
relationship_order <- c(
  "parent_child", "full_siblings", "half_siblings",
  "cousins", "second_cousins", "unrelated"
)
relationship_labels <- c(
  "parent_child"   = "Parent-Child",
  "full_siblings"  = "Full Siblings",
  "half_siblings"  = "Half Siblings",
  "cousins"        = "Cousins",
  "second_cousins" = "Second Cousins",
  "unrelated"      = "Unrelated"
)

loci_set_order <- c(
  "core_13", "identifiler_15", "expanded_20", "supplementary", "autosomal_29"
)
loci_set_labels <- c(
  "core_13"        = "Core 13",
  "identifiler_15" = "Identifiler 15",
  "expanded_20"    = "Expanded 20",
  "supplementary"  = "Supplementary",
  "autosomal_29"   = "Autosomal 29"
)

population_labels <- c(
  "AfAm"     = "African American",
  "Asian"    = "Asian",
  "Cauc"     = "Caucasian",
  "Hispanic" = "Hispanic",
  "all"      = "All Populations"
)

hypothesis_labels <- c(
  "parent_child"  = "Parent-Child Test",
  "full_siblings" = "Full Siblings Test"
)

# Population factor order — "all" placed last so it appears as the rightmost
# facet column, visually separated from the four named populations
population_order <- c("AfAm", "Asian", "Cauc", "Hispanic", "all")

# Color palettes — kept consistent with plots_matched_publication.R throughout
relationship_colors <- c(
  "Parent-Child"   = "#D55E00",   # Vermillion
  "Full Siblings"  = "#E69F00",   # Orange
  "Half Siblings"  = "#56B4E9",   # Sky blue
  "Cousins"        = "#009E73",   # Bluish green
  "Second Cousins" = "#CC79A7",   # Reddish purple
  "Unrelated"      = "#999999"    # Gray
)

population_colors <- c(
  "AfAm"     = "#0072B2",   # Deep blue
  "Asian"    = "#009E73",   # Bluish green
  "Cauc"     = "#56B4E9",   # Sky blue
  "Hispanic" = "#CC79A7",   # Reddish purple
  "all"      = "#999999"    # Gray (combined)
)


# Classification outcome colors
classification_colors <- c(
  "True Positive" = "#2ca02c",   # Green
  "Related FP"    = "#ff7f0e",   # Orange
  "Unrelated FP"  = "#d62728"    # Red
)

# Shared publication theme applied to both figures for visual consistency
theme_publication <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title      = element_text(hjust = 0.5, size = base_size + 3),
      plot.subtitle   = element_text(hjust = 0.5, size = base_size + 1, color = "grey40"),
      plot.caption    = element_text(hjust = 0, size = base_size - 2, color = "grey50"),
      axis.title      = element_text(size = base_size + 1),
      axis.text.x     = element_text(angle = 45, hjust = 1, size = base_size - 1),
      axis.text.y     = element_text(size = base_size - 1),
      strip.text      = element_text(size = base_size, face = "plain"),
      legend.title    = element_text(size = base_size),
      legend.text     = element_text(size = base_size - 1),
      legend.position = "bottom",
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.border       = element_rect(color = "grey85", fill = NA, linewidth = 0.5)
    )
}


# =============================================================================
# SECTION 2: LOAD DATA
# =============================================================================

log_message("Loading proportions data...")

proportions_file <- file.path(input_dir, "proportions_with_classification.csv")

if (!file.exists(proportions_file)) {
  stop(sprintf(
    "Intermediate file not found: %s\nRun prepare_combined_lr_intermediates.R first.",
    proportions_file
  ))
}

proportions_all <- fread(proportions_file) %>%
  mutate(
    known_relationship  = factor(known_relationship,  levels = relationship_order),
    tested_relationship = factor(tested_relationship, levels = relationship_order),
    loci_set            = factor(loci_set,            levels = loci_set_order),
    population          = factor(population,          levels = population_order),
    classification      = factor(classification,
                                 levels = c("True Positive", "Related FP", "Unrelated FP"))
  )

log_message(paste("Loaded", nrow(proportions_all), "rows."))


# =============================================================================
# SECTION 3: FIGURE 1 — MAIN TEXT
# Classification performance at Core 13 vs Expanded 20 loci, 0.01% FPR only
# =============================================================================
# Single-threshold view using only the most stringent, forensically relevant
# cutoff. Both loci panels are shown in ONE faceted plot so that the two loci
# conditions for a given hypothesis sit as ADJACENT rows — making the Core 13
# vs Expanded 20 comparison a single-row glance rather than a cross-panel jump.
#
# Facet layout: rows = tested_relationship (outer) + loci_set (inner),
#               cols = known_relationship.
# This yields 4 row facets (2 hypotheses x 2 loci sets) and 6 column facets,
# with a single shared set of column headers and one x-axis.
#
# Y-axis is FIXED at 0-100% across all facets (both loci sets reach ~100% at the
# True Positive / Related FP cells), so bar HEIGHTS are directly comparable
# everywhere — no free_y.
#
# The previous stacked A/B version is preserved, commented out, at the end of
# this section in case panel-letter referencing is preferred.
# =============================================================================

log_message("Building Figure 1: classification summary at 0.01% FPR, Core 13 vs Expanded 20...")

fig1_data <- proportions_all %>%
  filter(
    loci_set %in% c("core_13", "expanded_20"),
    tested_relationship %in% c("parent_child", "full_siblings")
  )

# prop_LR_gt_1000 = proportion of pairs exceeding the 0.01% FPR cutoff
figure1 <- ggplot(
  fig1_data,
  aes(
    x    = classification,
    y    = prop_LR_gt_1000,
    fill = population
  )
) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8) +
  # Rows: hypothesis (outer) + loci set (inner, so the two loci conditions are
  # adjacent rows); Cols: known relationship. Fixed y for direct comparability.
  facet_grid(
    tested_relationship + loci_set ~ known_relationship,
    labeller = labeller(
      tested_relationship = as_labeller(hypothesis_labels),
      loci_set            = as_labeller(loci_set_labels),
      known_relationship  = as_labeller(relationship_labels)
    )
  ) +
  scale_fill_manual(values = population_colors, name = "Population",
                    labels = population_labels) +
  scale_y_continuous(
    labels = percent_format(accuracy = 0.001),
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25)
  ) +
  labs(
    title    = "LR Classification Performance: Core 13 vs Expanded 20 Loci",
    subtitle = "Threshold: 0.01% false positive rate among unrelated pairs",
    x        = "Classification",
    y        = "Proportion of pairs with LR > threshold"
  ) +
  theme_publication(base_size = 14)

ggsave(
  filename = file.path(output_dir, "cutoff_classification_0.01fpr_core13_expanded20.png"),
  plot     = figure1,
  width    = 18,
  height   = 16,
  dpi      = 300,
  bg       = "white"
)

log_message("Figure 1 saved.")

# ---------------------------------------------------------------------------
# ALTERNATIVE (previous) LAYOUT: stacked A/B sub-panels, one per loci set.
# Uncomment this block — and comment out the faceted `figure1` above — to
# switch back to panel-letter (A/B) referencing with per-panel free y-axes.
# ---------------------------------------------------------------------------
# make_fig1_panel <- function(data, loci, panel_title) {
#   panel_data <- data %>% filter(loci_set == loci)
#   ggplot(
#     panel_data,
#     aes(x = classification, y = prop_LR_gt_1000, fill = population)
#   ) +
#     geom_col(position = position_dodge(width = 0.85), width = 0.8) +
#     facet_grid(
#       tested_relationship ~ known_relationship,
#       scales   = "free_y",
#       labeller = labeller(
#         tested_relationship = as_labeller(hypothesis_labels),
#         known_relationship  = as_labeller(relationship_labels)
#       )
#     ) +
#     scale_fill_manual(values = population_colors, name = "Population",
#                       labels = population_labels) +
#     scale_y_continuous(labels = percent_format(accuracy = 0.001)) +
#     labs(title = panel_title, x = "Classification",
#          y = "Proportion of pairs with LR > threshold") +
#     theme_publication(base_size = 14)
# }
#
# fig1a <- make_fig1_panel(fig1_data, "core_13",     "Core 13 Loci")
# fig1b <- make_fig1_panel(fig1_data, "expanded_20", "Expanded 20 Loci")
#
# figure1_ab <- (fig1a / fig1b) +
#   plot_annotation(
#     title      = "LR Classification Performance: Core 13 vs Expanded 20 Loci",
#     subtitle   = "Threshold: 0.01% false positive rate among unrelated pairs",
#     tag_levels = "A"
#   ) +
#   plot_layout(guides = "collect") &
#   theme(legend.position = "bottom")
#
# ggsave(
#   filename = file.path(output_dir, "cutoff_classification_0.01fpr_core13_expanded20_AB.png"),
#   plot     = figure1_ab,
#   width    = 18, height = 20, dpi = 300, bg = "white"
# )


# =============================================================================
# SECTION 4: FIGURE 2 — SUPPLEMENT
# Classification performance at 29 loci across all three FPR thresholds
# =============================================================================
# Extends Figure 1 by pivoting all three FPR thresholds into a row facet so
# the reader can see how classification changes with threshold stringency.
# 6 row facets total (2 hypotheses x 3 thresholds).
# Fixed (LR > 1) excluded — FPR-derived thresholds are the forensically
# relevant comparisons.
# =============================================================================

log_message("Building Figure 2: classification summary across all FPR thresholds, 29 loci...")

fig2_data <- proportions_all %>%
  filter(
    loci_set == "autosomal_29",
    tested_relationship %in% c("parent_child", "full_siblings")
  ) %>%
  pivot_longer(
    cols      = c(prop_LR_gt_10, prop_LR_gt_100, prop_LR_gt_1000),
    names_to  = "threshold",
    values_to = "proportion"
  ) %>%
  mutate(
    threshold = factor(
      threshold,
      levels = c("prop_LR_gt_10", "prop_LR_gt_100", "prop_LR_gt_1000"),
      labels = c("1% FPR", "0.1% FPR", "0.01% FPR")
    )
  )

figure2 <- ggplot(
  fig2_data,
  aes(
    x    = known_relationship,
    y    = proportion,
    fill = classification
  )
) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8) +
  # Row facets: tested hypothesis x threshold (6 rows total)
  # Column facets: population (5 columns including "all" as rightmost)
  facet_grid(
    tested_relationship + threshold ~ population,
    labeller = labeller(
      tested_relationship = as_labeller(hypothesis_labels),
      population          = as_labeller(population_labels)
    )
  ) +
  scale_fill_manual(values = classification_colors, name = "Classification") +
  scale_x_discrete(labels = relationship_labels) +
  scale_y_continuous(
    labels = percent_format(accuracy = 0.001),
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25)
  ) +
  labs(
    title    = "LR Classification Performance at 29 Autosomal Loci",
    subtitle = "Using population-matched allele frequencies across three false positive rate thresholds",
    x        = "True Relationship",
    y        = "Proportion of pairs with LR > threshold"
  ) +
  theme_publication(base_size = 14)

ggsave(
  filename = file.path(output_dir, "cutoff_classification_fpr_29loci.png"),
  plot     = figure2,
  width    = 18,
  height   = 16,
  dpi      = 300,
  bg       = "white"
)

log_message("Figure 2 saved.")


# =============================================================================
# SECTION 5: FIGURE 3 — SUPPLEMENT
# Heatmap of false positive rates across all loci panels and FPR cutoffs
# =============================================================================
# Heatmap encodes false positive rate as color intensity across all loci panels,
# stratified by population and tested hypothesis. Restricted to false positive
# rows only so the color scale is entirely dedicated to showing FP severity.
# All four thresholds (including Fixed LR > 1) shown for completeness.
# =============================================================================

log_message("Building Figure 3: supplement heatmap of false positive rates...")

fig3_data <- proportions_all %>%
  filter(
    tested_relationship %in% c("parent_child", "full_siblings"),
    classification != "True Positive",   # focus exclusively on false positive rates
    population != "all",                  # exclude pooled category before averaging
    tested_population == population    # add this line
  ) %>%
  # Pivot all four threshold columns so cutoff becomes a faceting variable
  pivot_longer(
    cols      = c(prop_LR_gt_1, prop_LR_gt_10, prop_LR_gt_100, prop_LR_gt_1000),
    names_to  = "threshold",
    values_to = "proportion"
  ) %>%
  mutate(
    threshold = factor(
      threshold,
      levels = c("prop_LR_gt_1", "prop_LR_gt_10", "prop_LR_gt_100", "prop_LR_gt_1000"),
      labels = c("Fixed (LR > 1)", "1% FPR", "0.1% FPR", "0.01% FPR")
    )
   ) #%>%
  # # Average proportion across populations for each cell
  # group_by(known_relationship, loci_set, tested_relationship, threshold, classification) %>%
  # summarise(mean_proportion = mean(proportion, na.rm = TRUE), .groups = "drop")

make_heatmap <- function(data, hypothesis, panel_title) {
  panel_data <- data %>% filter(tested_relationship == hypothesis)
  # Relative threshold for white text: top 45% of this panel's own scale
  color_threshold <- max(panel_data$proportion, na.rm = TRUE) * 0.55

  ggplot(panel_data, aes(x = loci_set, y = known_relationship, fill = proportion)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(
      aes(
        label = percent(proportion, accuracy = 0.001),
        color = proportion > color_threshold
      ),
      size = 2.8
    ) +
    scale_color_manual(values = c("TRUE" = "white", "FALSE" = "grey20"), guide = "none") +
    facet_grid(
      population ~ threshold,
      labeller = labeller(
        population = as_labeller(population_labels)
      )
    ) +
    scale_fill_distiller(
      palette   = "YlOrRd",
      direction = 1,
      labels    = percent_format(accuracy = 0.001),
      name      = "False positive rate"
    ) +
    scale_x_discrete(labels = loci_set_labels) +
    scale_y_discrete(labels = relationship_labels) +
    labs(x = "Loci Panel", y = "True Relationship", title = panel_title) +
    theme_publication(base_size = 14) +
    theme(
      panel.grid.major  = element_blank(),
      axis.text.x       = element_text(angle = 30, hjust = 1),
      legend.text       = element_text(angle = 45, hjust = 1),
      legend.key.width  = unit(1.5, "cm"),
      legend.title.position = "top"
    )
}


fig3a <- make_heatmap(fig3_data, "parent_child", "Parent-Child Test")
fig3b <- make_heatmap(fig3_data, "full_siblings", "Full Siblings Test")

figure3 <- (fig3a / fig3b) +
  plot_annotation(
    title      = "False Positive Rates Across Loci Panels and LR Thresholds",
    tag_levels = "A"
  ) &
  theme(legend.position = "bottom")

ggsave(
  filename = file.path(output_dir, "cutoff_supp_heatmap_fp_rates.png"),
  plot     = figure3,
  width    = 22,
  height   = 20,
  dpi      = 300,
  bg       = "white"
)

log_message("Figure 3 saved.")
log_message("All manuscript figures complete.")