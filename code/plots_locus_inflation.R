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
#   Fig 2:  Combined A/B figure:
#     A)    All mismatch directions — heatmap of per-locus median log10 ratio
#           for full siblings across all true x tested population combinations.
#     B)    Observed heterozygosity for primary driver loci by population.
#
#   In all figures loci are arranged along their axis by nested CODIS set
#   (Core 13, then the additions that make Identifiler 15, Expanded 20,
#   Supplementary, then Autosomal-29-only loci), with separator lines and
#   group brackets/strips marking each set. Within a set, loci are ordered by
#   AfAm->Asian full-sibling inflation (see WITHIN_GROUP_ORDER).
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
                 "Supplementary", "Autosomal 29")
PANEL_COLORS <- c(
  "Core 13"          = "#E69F00",
  "Identifiler 15"   = "#56B4E9",
  "Expanded 20"      = "#009E73",
  "Supplementary"    = "#CC79A7",
  "Autosomal 29"= "#D55E00"
)

# Short panel labels for compact group brackets on the figures
PANEL_SHORT <- c(
  "Core 13"           = "Core 13",
  "Identifiler 15"    = "Ident. 15",
  "Expanded 20"       = "Exp. 20",
  "Supplementary"     = "Suppl.",
  "Autosomal 29" = "Auto. 29"
)

# How to order loci *within* each CODIS group along the axis:
#   "inflation"  -> descending AfAm->Asian full-sib inflation (keeps driver rank
#                   visible inside each set; this is the default)
#   "alpha"      -> alphabetical by locus name
WITHIN_GROUP_ORDER <- "inflation"

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

# -----------------------------------------------------------------------------
# axis_group_geometry()
# Given the loci as they appear along an axis (position 1 = first level) and a
# named vector mapping locus -> CODIS panel, return one row per contiguous
# CODIS group with its start / end / midpoint position. Used to draw the
# separator lines and the group brackets/labels.
# -----------------------------------------------------------------------------
axis_group_geometry <- function(levels_in_order, panel_map) {
  panels <- as.character(panel_map[levels_in_order])
  r      <- rle(panels)
  ends   <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L
  tibble::tibble(
    panel = r$values,
    start = starts,
    end   = ends,
    mid   = (starts + ends) / 2
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
      TRUE                ~ "Autosomal 29"
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
               panel = factor("Autosomal 29", levels = PANEL_ORDER))
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

# -----------------------------------------------------------------------------
# CANONICAL LOCUS ORDER — nested CODIS sets
#
# Loci are grouped along the axis following the nested CODIS build:
#   Core 13  ->  (+) Identifiler 15  ->  (+) Expanded 20  ->
#   (+) Supplementary  ->  Autosomal 29
#
# Within each CODIS group loci are ordered by WITHIN_GROUP_ORDER:
#   "inflation" (default) uses the canonical AfAm->Asian full-sib inflation so
#   the driver ranking is still visible inside each set; loci without a value
#   fall to the end of their group. "alpha" sorts by locus name.
#
# codis_order_top_to_bottom is the intended top-to-bottom reading order, with
# Core 13 at the TOP. Because ggplot draws the first factor level at the bottom
# of a y axis, locus_order (used by Figs 1 and 2A) is its reverse.
# -----------------------------------------------------------------------------

# Per-locus value used to order loci within a CODIS group
within_group_rank <- inflation %>%
  filter(known_relationship == "Full Siblings",
         population == "AfAm", tested_population == "Asian") %>%
  select(locus, rank_val = median_log10_ratio)

codis_order_top_to_bottom <- panel_lookup %>%
  distinct(locus, panel) %>%
  filter(locus %in% all_loci) %>%                        # only loci present in data
  left_join(within_group_rank, by = "locus") %>%
  mutate(
    panel    = factor(panel, levels = PANEL_ORDER),
    rank_val = if (WITHIN_GROUP_ORDER == "alpha") 0 else replace_na(rank_val, -Inf)
  ) %>%
  arrange(panel, desc(rank_val), locus) %>%
  pull(locus)

# locus -> panel lookup as a named vector, for separators / brackets
locus_panel <- panel_lookup %>%
  distinct(locus, panel) %>%
  filter(locus %in% all_loci) %>%
  deframe()

# Factor levels for the y axes (Core 13 ends up at the top)
locus_order <- rev(codis_order_top_to_bottom)

fig1_data <- fig1_data %>%
  mutate(locus_ordered = factor(locus, levels = locus_order),
         ypos          = as.integer(locus_ordered))  # numeric y for brackets

# CODIS group geometry along the y axis (positions match locus_order)
fig1_groups     <- axis_group_geometry(locus_order, locus_panel)
fig1_boundaries <- head(fig1_groups$end, -1) + 0.5

# Right-side group brackets are drawn only in the rightmost relationship facet
rightmost_rel <- tail(levels(droplevels(fig1_data$known_relationship)), 1)
x_rng  <- range(c(fig1_data$ci_lower, fig1_data$ci_upper,
                  fig1_data$median_log10_ratio), na.rm = TRUE)
x_span <- diff(x_rng); if (!is.finite(x_span) || x_span == 0) x_span <- 1
x_br   <- x_rng[2] + 0.03 * x_span   # vertical bracket bar
x_tick <- 0.012 * x_span             # bracket tick length
x_txt  <- x_rng[2] + 0.055 * x_span  # label position

fig1_brackets <- tidyr::crossing(fig1_groups,
                                  direction = unique(fig1_data$direction)) %>%
  mutate(known_relationship = factor(rightmost_rel,
                                     levels = levels(fig1_data$known_relationship)))

fig1 <- ggplot(fig1_data,
               aes(x = median_log10_ratio, y = ypos, color = panel)) +

  # Zero reference line
  geom_vline(xintercept = 0, color = "gray60", linewidth = 0.4, linetype = "dashed") +

  # Driver threshold line
  geom_vline(xintercept = DRIVER_THRESHOLD, color = "gray40", linewidth = 0.4,
             linetype = "dotted") +

  # CODIS group separators
  geom_hline(yintercept = fig1_boundaries, color = "gray75", linewidth = 0.4) +

  # CI segment
  geom_segment(aes(x = ci_lower, xend = ci_upper,
                   y = ypos, yend = ypos),
               alpha = 0.3, linewidth = 0.6) +

  # Point
  geom_point(size = 2.5) +

  # CODIS group brackets + labels (rightmost facet only, drawn in the margin)
  geom_segment(data = fig1_brackets, inherit.aes = FALSE,
               aes(y = start - 0.4, yend = end + 0.4),
               x = x_br, xend = x_br, color = "gray30", linewidth = 0.5) +
  geom_segment(data = fig1_brackets, inherit.aes = FALSE,
               aes(y = start - 0.4, yend = start - 0.4),
               x = x_br - x_tick, xend = x_br, color = "gray30", linewidth = 0.5) +
  geom_segment(data = fig1_brackets, inherit.aes = FALSE,
               aes(y = end + 0.4, yend = end + 0.4),
               x = x_br - x_tick, xend = x_br, color = "gray30", linewidth = 0.5) +
  geom_text(data = fig1_brackets, inherit.aes = FALSE,
            aes(y = mid, label = panel),
            x = x_txt, hjust = 0, size = 3.0, color = "gray20") +

  facet_grid(direction ~ known_relationship) +

  scale_color_manual(values = PANEL_COLORS, name = "CODIS Panel",
                     drop = FALSE) +

  scale_y_continuous(breaks = seq_along(locus_order), labels = locus_order,
                     expand = expansion(add = 0.6)) +

  coord_cartesian(xlim = c(x_rng[1] - 0.02 * x_span, x_rng[2] + 0.02 * x_span),
                  clip = "off") +

  labs(
    title    = "Per-Locus LR Inflation Under African American \u2194 Asian Frequency Mismatch",
    subtitle = paste0("Median log\u2081\u2080(LR_wrong / LR_correct) per locus | True positives only\n",
                      "Loci grouped by nested CODIS set (Core 13 at top) | ",
                      "dotted line = driver threshold (", DRIVER_THRESHOLD, ")"),
    x        = "Median log\u2081\u2080(LR_wrong / LR_correct)",
    y        = NULL
  ) +

  theme_publication() +
  theme(
    legend.position = "right",
    axis.text.y     = element_text(size = rel(0.75)),
    plot.margin     = margin(10, 115, 10, 10)   # room for right-side brackets
  )

fig1_file <- "locus_inflation_ranked.png"
ggsave(file.path(output_dir, fig1_file),
       plot = fig1, width = 12, height = 10, dpi = 300, bg = "white")
log_message(paste("Figure 1 saved:", fig1_file))

# =============================================================================
# IDENTIFY PRIMARY DRIVER LOCI — direction-specific
# Each direction has its own ordered loci set for Panel B.
# Union used for the driver loci table.
# =============================================================================

# AfAm->Asian: loci above threshold, ordered by that direction's inflation
driver_loci_afam_asian <- inflation %>%
  filter(known_relationship == "Full Siblings",
         population == "AfAm", tested_population == "Asian",
         median_log10_ratio >= DRIVER_THRESHOLD) %>%
  arrange(desc(median_log10_ratio)) %>%
  pull(locus)

# Asian->AfAm: loci above threshold, ordered by that direction's inflation
driver_loci_asian_afam <- inflation %>%
  filter(known_relationship == "Full Siblings",
         population == "Asian", tested_population == "AfAm",
         median_log10_ratio >= DRIVER_THRESHOLD) %>%
  arrange(desc(median_log10_ratio)) %>%
  pull(locus)

# Union for the driver loci table
driver_loci <- union(driver_loci_afam_asian, driver_loci_asian_afam)

log_message(sprintf(
  "Driver loci: %d in AfAm->Asian, %d in Asian->AfAm, %d in overlap, %d total.",
  length(driver_loci_afam_asian), length(driver_loci_asian_afam),
  length(intersect(driver_loci_afam_asian, driver_loci_asian_afam)),
  length(driver_loci)
))

# =============================================================================
# FIGURE 2A: HEATMAP — all mismatch directions, full siblings
# =============================================================================

log_message("Generating Figure 2A: All-direction mismatch heatmap...")

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
    # Same locus order as Fig 1; the nested CODIS set becomes a facet strip
    locus_ordered = factor(locus, levels = locus_order),
    panel         = factor(panel, levels = PANEL_ORDER)
  )

fig2 <- ggplot(fig2_data,
               aes(x = tested_pop, y = locus_ordered,
                   fill = median_log10_ratio)) +

  geom_tile(color = "white", linewidth = 0.4) +

  # Rows = nested CODIS set (labelled strip acts as the group bracket);
  # free/space y keeps tiles uniform and sizes each set by its locus count.
  facet_grid(panel ~ true_pop, scales = "free_y", space = "free_y") +

  scale_fill_gradientn(
    colors = c("white", "#c6dbef", "#6baed6", "#2166ac", "#08306b"),
    limits = c(0, NA),
    oob    = scales::squish,
    name   = "Median\nlog\u2081\u2080(LR_wrong /\nLR_correct)"
  ) +

  labs(
    title    = "Per-Locus LR Inflation Across All Population Mismatch Directions",
    subtitle = "Full siblings | Loci grouped by nested CODIS set | each cell = one locus x tested frequency database",
    x        = "Tested Population Frequencies",
    y        = NULL
  ) +

  theme_publication() +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = rel(0.8)),
    axis.text.y     = element_text(size  = rel(0.7)),
    panel.grid      = element_blank(),
    panel.spacing.y = unit(4, "pt"),
    strip.text.y    = element_text(angle = 0, size = rel(0.75)),
    legend.position = "right"
  )

fig2_file <- "locus_inflation_heatmap_heterozygosity.png"  # combined output

# =============================================================================
# FIGURE 2B: OBSERVED HETEROZYGOSITY — all driver loci, single panel
#
# All driver loci shown. Color = population. Shape = which mismatch
# direction(s) the locus drives: AfAm->Asian only, Asian->AfAm only, or Both.
# Loci ordered left to right by nested CODIS set (Core 13 first).
# =============================================================================

log_message("Generating Figure 2B: Observed heterozygosity at driver loci...")

# Order driver loci by nested CODIS set (Core 13 leftmost), consistent with
# the canonical order used in Figs 1 and 2A.
driver_loci_ordered <- codis_order_top_to_bottom[
  codis_order_top_to_bottom %in% driver_loci]

# Assign direction category per locus
direction_category <- tibble(locus = driver_loci) %>%
  mutate(
    driver_direction = case_when(
      locus %in% driver_loci_afam_asian & locus %in% driver_loci_asian_afam ~ "Both directions",
      locus %in% driver_loci_afam_asian                                      ~ "AfAm\u2192Asian only",
      TRUE                                                                    ~ "Asian\u2192AfAm only"
    ),
    driver_direction = factor(driver_direction,
                              levels = c("Both directions",
                                         "AfAm\u2192Asian only",
                                         "Asian\u2192AfAm only"))
  )

# Long format for plotting
heteroz_long <- heteroz %>%
  filter(locus %in% driver_loci_ordered) %>%
  select(locus, AfAm, Asian, Cauc, Hispanic) %>%
  pivot_longer(cols      = c(AfAm, Asian, Cauc, Hispanic),
               names_to  = "population",
               values_to = "obs_heterozygosity") %>%
  left_join(direction_category, by = "locus") %>%
  mutate(
    population    = factor(population, levels = names(POPULATION_LABELS),
                           labels = POPULATION_LABELS),
    locus_ordered = factor(locus, levels = driver_loci_ordered),
    xpos          = as.integer(locus_ordered)   # numeric x for brackets
  )

# CODIS group geometry along the x axis (positions match driver_loci_ordered)
fig3_groups     <- axis_group_geometry(driver_loci_ordered, locus_panel)
fig3_boundaries <- head(fig3_groups$end, -1) + 0.5
fig3_brackets   <- fig3_groups %>%
  mutate(label = PANEL_SHORT[panel])

y_br  <- 1.035   # bracket bar height (inside an expanded header band, not the margin)
y_tk  <- 0.012   # bracket tick depth
y_txt <- 1.05    # label baseline

fig3 <- ggplot(heteroz_long,
               aes(x     = xpos,
                   y     = obs_heterozygosity,
                   color = population,
                   shape = driver_direction,
                   group = population)) +

  # CODIS group separators
  geom_vline(xintercept = fig3_boundaries, color = "gray80", linewidth = 0.4) +

  geom_line(alpha = 0.5, linewidth = 0.6) +
  geom_point(size = 3) +

  # CODIS group brackets + labels above the panel
  geom_segment(data = fig3_brackets, inherit.aes = FALSE,
               aes(x = start - 0.4, xend = end + 0.4), y = y_br, yend = y_br,
               color = "gray30", linewidth = 0.5) +
  geom_segment(data = fig3_brackets, inherit.aes = FALSE,
               aes(x = start - 0.4, xend = start - 0.4),
               y = y_br, yend = y_br - y_tk, color = "gray30", linewidth = 0.5) +
  geom_segment(data = fig3_brackets, inherit.aes = FALSE,
               aes(x = end + 0.4, xend = end + 0.4),
               y = y_br, yend = y_br - y_tk, color = "gray30", linewidth = 0.5) +
  geom_text(data = fig3_brackets, inherit.aes = FALSE,
            aes(x = mid, label = label), y = y_txt,
            vjust = 0, size = 3.0, color = "gray20") +

  scale_color_manual(values = setNames(POPULATION_COLORS,
                                       POPULATION_LABELS[names(POPULATION_COLORS)]),
                     name = "Population") +
  scale_shape_manual(values = c("Both directions"  = 16,   # filled circle
                                "AfAm\u2192Asian only" = 17,   # filled triangle
                                "Asian\u2192AfAm only" = 15),  # filled square
                     name = "Driver Direction") +
  scale_x_continuous(breaks = seq_along(driver_loci_ordered),
                     labels = driver_loci_ordered,
                     expand = expansion(add = 0.6)) +
  scale_y_continuous(breaks = seq(0.6, 1.0, 0.1),
                     labels = scales::percent_format(accuracy = 1)) +

  # Expand the upper limit to open a header band above 100% for the CODIS
  # brackets, so they sit inside the panel and clear of the subtitle.
  coord_cartesian(ylim = c(0.55, 1.09), clip = "off") +

  labs(
    title    = NULL,
    subtitle = paste0("Primary driver loci (median log\u2081\u2080(LR_wrong / LR_correct) \u2265 ",
                      DRIVER_THRESHOLD, " in AfAm \u2194 Asian full siblings)\n",
                      "Ordered by nested CODIS set"),
    x        = "Locus",
    y        = "Observed Heterozygosity"
  ) +

  theme_publication() +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.margin     = margin(10, 10, 10, 10)
  )

# =============================================================================
# COMBINE 2A + 2B WITH PATCHWORK AND SAVE
# =============================================================================

log_message("Combining panels and saving Figure 2 (A/B)...")

fig2_combined <- fig2 / fig3 +
  plot_layout(heights = c(3, 2)) +
  plot_annotation(
    tag_levels = "A",
    theme      = theme(plot.tag = element_text(size = 16, face = "bold"))
  )

ggsave(file.path(output_dir, fig2_file),
       plot = fig2_combined, width = 14, height = 18, dpi = 300, bg = "white")
log_message(paste("Figure 2 (A/B) saved:", fig2_file))

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
    driver_AfAm_Asian = locus %in% driver_loci_afam_asian,
    driver_Asian_AfAm = locus %in% driver_loci_asian_afam,
    is_primary_driver = driver_AfAm_Asian | driver_Asian_AfAm,
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
for (f in c(fig1_file, fig2_file, "primary_driver_loci_table.csv")) {
  cat(sprintf("  %s\n", file.path(output_dir, f)))
}
cat("=============================================================================\n")

log_message("Done.")
