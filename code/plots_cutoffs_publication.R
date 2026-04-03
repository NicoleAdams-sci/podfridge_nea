#!/usr/bin/env Rscript

# =============================================================================
# Manuscript Figures: LR Classification Performance
# =============================================================================
# Produces three figures:
#
#   FIGURE 1 (Main Text):
#     Classification summary at 29 autosomal loci using the 0.01% FPR threshold
#     only. Clean, single-threshold view for the main text argument.
#     Output: cutoff_classification_0.01fpr_29loci.png
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
# Input:  proportions_with_classification.csv (output of plots_proportion_exceeding_cutoffs.R)
#
# Usage:
#   Rscript plots_cutoffs_publication.R <input_dir> [output_dir]
# =============================================================================

suppressMessages(suppressWarnings({
  library(tidyverse)
  library(scales)
  library(patchwork)
}))

log_message <- function(msg) cat(paste0("[", Sys.time(), "] ", msg, "\n"))

# --- Argument parsing ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript plots_manuscript_figures.R <input_dir> [output_dir]")

input_dir  <- args[1]
output_dir <- if (length(args) >= 2) file.path("output", args[2]) else file.path("output", input_dir)
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

# Color palettes — kept consistent with plots_matched.R throughout
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
# SECTION 2: LOAD AND PREPARE DATA
# =============================================================================

log_message("Loading proportions data...")

proportions_all <- read_csv(
  file.path(input_dir, "plots_exceeding_cutoffs/proportions_with_classification.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    known_relationship  = factor(known_relationship,  levels = relationship_order),
    tested_relationship = factor(tested_relationship, levels = relationship_order),
    loci_set            = factor(loci_set,            levels = loci_set_order),
    population          = factor(population,          levels = population_order),
    classification      = factor(
      classification,
      levels = c("True Positive", "Related FP", "Unrelated FP")
    )
  )

log_message(paste("Loaded", nrow(proportions_all), "rows."))


# =============================================================================
# SECTION 3: FIGURE 1 — MAIN TEXT
# Classification performance at 29 loci, 0.01% FPR threshold only
# =============================================================================
# Single-threshold view using only the most stringent, forensically relevant
# cutoff. Two row facets (one per hypothesis), five column facets (four
# populations + "all" rightmost). Kept clean for the main text argument.
# =============================================================================

log_message("Building Figure 1: classification summary at 0.01% FPR only, 29 loci...")

fig1_data <- proportions_all %>%
  filter(
    loci_set == "autosomal_29",
    tested_relationship %in% c("parent_child", "full_siblings")
  )

figure1 <- ggplot(
  fig1_data,
  aes(
    x    = classification,
    y    = prop_LR_gt_1000,     # prop_LR_gt_1000 = proportion exceeding 0.01% FPR cutoff
    fill = population
  )
) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8) +
  facet_grid(
    tested_relationship ~ known_relationship,
    scales = "free_y",
    labeller = labeller(
      tested_relationship = as_labeller(hypothesis_labels),
      known_relationship          = as_labeller(relationship_labels)
    )
  ) +
  scale_fill_manual(values = population_colors, name = "Population",
                    labels = population_labels) +
  #scale_x_discrete(labels = relationship_labels) +
  scale_y_continuous(
    labels = percent_format(accuracy = 0.001)
  ) +
  labs(
    title    = "LR Classification Performance at 29 Autosomal Loci",
    subtitle = "Threshold: 0.01% false positive rate among unrelated pairs",
    x        = "Classification",
    y        = "Proportion of pairs with LR > threshold"
  ) +
  theme_publication(base_size = 14)

ggsave(
  filename = file.path(output_dir, "cutoff_classification_0.01fpr_29loci.png"),
  plot     = figure1,
  width    = 18,
  height   = 10,
  dpi      = 300,
  bg       = "white"
)

log_message("Figure 1 saved.")


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
  ggplot(data %>% filter(tested_relationship == hypothesis),
         aes(x = loci_set, y = known_relationship, fill = proportion)
  ) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(
      aes(
        label = percent(proportion, accuracy = 0.001),
        color = proportion > 0.55
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
      limits    = c(0, 1),
      labels    = percent_format(accuracy = 0.001),
      name      = "False positive rate"
    ) +
    scale_x_discrete(labels = loci_set_labels) +
    scale_y_discrete(labels = relationship_labels) +
    labs(x = "Loci Panel", y = "True Relationship", title = panel_title) +
    theme_publication(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      axis.text.x      = element_text(angle = 30, hjust = 1)
    )
}


fig3a <- make_heatmap(fig3_data, "parent_child", "Parent-Child Test")
fig3b <- make_heatmap(fig3_data, "full_siblings", "Full Siblings Test")

figure3 <- (fig3a / fig3b) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title      = "False Positive Rates Across Loci Panels and LR Thresholds",
    tag_levels = "A"
  ) &
  theme(legend.position = "bottom")   # the & applies to all panels including collected legend

ggsave(
  filename = file.path(output_dir, "cutoff_supp_heatmap_fp_rates.png"),
  plot     = figure3,
  width    = 16,
  height   = 20,
  dpi      = 300,
  bg       = "white"
)

log_message("Figure 3 saved.")
log_message("All manuscript figures complete.")