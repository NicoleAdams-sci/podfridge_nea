#!/usr/bin/env Rscript
# ==============================================================================
# Publication-Ready Plots for Population & Relationship Matched LR Analysis
# ==============================================================================
# Creates:
#   1. Main text figure: Violin plots by relationship type (all populations)
#   2. Supplement figure: Box plots by population (showing consistency)
#   3. Statistical summary table (CSV) for results section
#
# Statistical tests:
#   - Kruskal-Wallis + epsilon-squared effect size: LR ~ relationship type
#     (per loci set, "all" population) 
#   - Spearman correlation: median log10(LR) ~ loci count (per relationship)
#   - Kruskal-Wallis + epsilon-squared: population effect at Autosomal 29
#     (per relationship type)
#
# NOTE: With large simulation sample sizes, p-values will be near-zero for
# most tests. Effect sizes (epsilon-squared) are the informative statistic.
#
# Date: 2026-03-06
# ==============================================================================

# Load Required Libraries ----
suppressMessages(suppressWarnings({
  library(tidyverse)
  library(data.table)
  library(scales)
  library(patchwork)
  library(rstatix)   # For kruskal_effsize (epsilon-squared)
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
combined_lrs_match <- fread(raw_lrs_path) %>%
  as_tibble() %>%
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
    axis.title.y       = element_text(size = 14, face = "bold"),
    strip.text         = element_text(size = 13, face = "bold"),
    strip.background   = element_rect(fill = "gray90", color = "gray60"),
    legend.position    = "bottom",
    legend.title       = element_text(size = 13, face = "bold"),
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
    axis.title       = element_text(size = 14, face = "bold"),
    strip.text       = element_text(size = 13, face = "bold"),
    strip.background = element_rect(fill = "gray90", color = "gray60"),
    legend.position  = "bottom",
    legend.title     = element_text(size = 13, face = "bold"),
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

# ==============================================================================
# STATISTICAL TESTS
# ==============================================================================
# Three targeted tests aligned with the matched results narrative:
#
#   TEST 1: Does log10(LR) differ by relationship type?
#           Kruskal-Wallis per loci set ("all" population only)
#           + epsilon-squared (ε²) effect size
#           ε² interpretation: small ≥ 0.01, medium ≥ 0.06, large ≥ 0.14
#
#   TEST 2: Does log10(LR) increase with more loci?
#           Spearman correlation of median log10(LR) vs. locus count,
#           per relationship type
#
#   TEST 3: Does population affect log10(LR) in the matched analysis?
#           Kruskal-Wallis per relationship type at Autosomal 29 only
#           + epsilon-squared effect size
#           (uses named populations only, excludes "all")
# ==============================================================================

log_message("Running statistical tests...")

# ------------------------------------------------------------------------------
# TEST 1: Relationship type effect per loci set
# ------------------------------------------------------------------------------
# Two-step approach to avoid cur_data() segfault with large datasets:
# Step 1a: Kruskal-Wallis test statistics via summarise
kw_rel_test <- main_data %>%
  group_by(loci_set) %>%
  summarise(
    kw_statistic = kruskal.test(log_LR ~ relationship)$statistic,
    kw_df        = kruskal.test(log_LR ~ relationship)$parameter,
    kw_p_value   = kruskal.test(log_LR ~ relationship)$p.value,
    n_pairs      = n(),
    .groups = "drop"
  )

# Step 1b: Epsilon-squared via group_modify (avoids cur_data() crash)
kw_rel_effsize <- main_data %>%
  group_by(loci_set) %>%
  group_modify(~ {
    es <- kruskal_effsize(data = .x, formula = log_LR ~ relationship, ci = FALSE)
    tibble(epsilon_sq = es$effsize, effect_size_magnitude = as.character(es$magnitude))
  }) %>%
  ungroup()

kw_relationship <- left_join(kw_rel_test, kw_rel_effsize, by = "loci_set") %>%
  mutate(
    test     = "Kruskal-Wallis",
    variable = "Relationship type",
    note     = "All populations combined (population = 'all')"
  )

cat("\n=== TEST 1: Relationship Type Effect (per loci set) ===\n")
print(kw_relationship %>%
  select(loci_set, kw_statistic, kw_df, kw_p_value, epsilon_sq, effect_size_magnitude))

# ------------------------------------------------------------------------------
# TEST 2: Loci count effect on median log10(LR) per relationship
# ------------------------------------------------------------------------------
# Compute median log10(LR) per relationship x loci set, then correlate with
# locus count. Uses "all" population data (same as main figure).
median_by_loci <- main_data %>%
  group_by(relationship, loci_set) %>%
  summarise(median_logLR = median(log_LR), .groups = "drop") %>%
  mutate(loci_count = loci_counts[as.character(loci_set)])

spearman_loci <- median_by_loci %>%
  group_by(relationship) %>%
  summarise(
    spearman_rho = cor(loci_count, median_logLR, method = "spearman"),
    spearman_p   = cor.test(loci_count, median_logLR,
                            method = "spearman", exact = FALSE)$p.value,
    n_loci_sets  = n(),
    .groups = "drop"
  ) %>%
  mutate(
    test     = "Spearman correlation",
    variable = "Loci count vs. median log10(LR)"
  )

cat("\n=== TEST 2: Loci Count Effect (Spearman, per relationship) ===\n")
print(spearman_loci %>% select(relationship, spearman_rho, spearman_p))

# ------------------------------------------------------------------------------
# TEST 3: Population effect at Autosomal 29
# ------------------------------------------------------------------------------
pop_data_29 <- combined_lrs_match %>%
  filter(loci_set == "Autosomal 29",
         population != "all")   # exclude pooled population

# Step 3a: Kruskal-Wallis test statistics
kw_pop_test <- pop_data_29 %>%
  group_by(relationship) %>%
  summarise(
    kw_statistic = kruskal.test(log_LR ~ population)$statistic,
    kw_df        = kruskal.test(log_LR ~ population)$parameter,
    kw_p_value   = kruskal.test(log_LR ~ population)$p.value,
    n_pairs      = n(),
    .groups = "drop"
  )

# Step 3b: Epsilon-squared via group_modify
kw_pop_effsize <- pop_data_29 %>%
  group_by(relationship) %>%
  group_modify(~ {
    es <- kruskal_effsize(data = .x, formula = log_LR ~ population, ci = FALSE)
    tibble(epsilon_sq = es$effsize, effect_size_magnitude = as.character(es$magnitude))
  }) %>%
  ungroup()

kw_population <- left_join(kw_pop_test, kw_pop_effsize, by = "relationship") %>%
  mutate(
    test     = "Kruskal-Wallis",
    variable = "Population (at Autosomal 29)",
    note     = "Named populations only (AfAm, Asian, Cauc, Hispanic)"
  )

cat("\n=== TEST 3: Population Effect at Autosomal 29 ===\n")
print(kw_population %>%
  select(relationship, kw_statistic, kw_df, kw_p_value, epsilon_sq, effect_size_magnitude))

# ------------------------------------------------------------------------------
# Compile and save statistical summary table
# ------------------------------------------------------------------------------
stats_table <- bind_rows(
  # Test 1
  kw_relationship %>%
    transmute(
      test,
      question    = "Does LR differ by relationship type?",
      grouping    = as.character(loci_set),
      statistic   = round(kw_statistic, 2),
      df          = kw_df,
      p_value     = formatC(kw_p_value, format = "e", digits = 2),
      epsilon_sq  = round(epsilon_sq, 3),
      magnitude   = effect_size_magnitude,
      note
    ),
  # Test 2
  spearman_loci %>%
    transmute(
      test,
      question    = "Does LR increase with more loci?",
      grouping    = as.character(relationship),
      statistic   = round(spearman_rho, 3),
      df          = NA_real_,
      p_value     = formatC(spearman_p, format = "e", digits = 2),
      epsilon_sq  = NA_real_,
      magnitude   = NA_character_,
      note        = "rho = Spearman rank correlation coefficient"
    ),
  # Test 3
  kw_population %>%
    transmute(
      test,
      question    = "Does population affect LR (matched analysis)?",
      grouping    = as.character(relationship),
      statistic   = round(kw_statistic, 2),
      df          = kw_df,
      p_value     = formatC(kw_p_value, format = "e", digits = 2),
      epsilon_sq  = round(epsilon_sq, 3),
      magnitude   = effect_size_magnitude,
      note        = "Autosomal 29 only; named populations (excludes 'all')"
    )
)

write_csv(stats_table,
          file.path(output_dir, "matched_statistical_tests.csv"))
log_message("Statistical test results saved to matched_statistical_tests.csv")

cat("\n=== FULL STATISTICAL SUMMARY TABLE ===\n")
print(stats_table, n = Inf)

# Descriptive Summary Statistics (unchanged from original)
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
