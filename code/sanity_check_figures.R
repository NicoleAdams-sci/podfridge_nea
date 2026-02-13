#!/usr/bin/env Rscript
# ==============================================================================
# sanity_check_figures.R
# Generate diagnostic figures from sanity check output
#
# Usage: Rscript code/sanity_check_figures.R [SANITY_CHECK_DIR]
# ==============================================================================

library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)

cat("=== Sanity Check Figures ===\n")

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1) args[1] else "output/sanity_check"

if (!dir.exists(input_dir)) {
  stop("Input directory not found: ", input_dir)
}

# Read the comprehensive results
results_file <- file.path(input_dir, "sanity_check_all_combinations.csv")
if (!file.exists(results_file)) {
  stop("Results file not found: ", results_file)
}

dt <- fread(results_file)
cat("Read", nrow(dt), "rows from", results_file, "\n")

# Define consistent ordering
rel_order <- c("parent_child", "full_siblings", "half_siblings",
               "cousins", "second_cousins", "unrelated")
rel_labels <- c("Parent-Child", "Full Siblings", "Half Siblings",
                "Cousins", "2nd Cousins", "Unrelated")
loci_order <- c("core_13", "identifiler_15", "expanded_20",
                "supplementary", "autosomal_29")
loci_labels <- c("Core 13", "Identifiler 15", "Expanded 20",
                 "Supplementary 23", "Autosomal 29")

dt[, known_relationship := factor(known_relationship, levels = rel_order, labels = rel_labels)]
dt[, tested_relationship := factor(tested_relationship,
     levels = c("parent_child", "full_siblings"),
     labels = c("Tested as Parent-Child", "Tested as Full Siblings"))]
dt[, loci_set := factor(loci_set, levels = loci_order, labels = loci_labels)]

# Keep all populations including "all"
dt_pop <- dt

# ==========================================================================
# Figure 1: Heatmap of FP rates for unrelated pairs (correct pop)
# ==========================================================================
cat("Generating Figure 1: Unrelated FP heatmap...\n")

fp_data <- dt_pop[classification == "Unrelated FP" & is_correct_pop == TRUE]

if (nrow(fp_data) > 0) {
  p1 <- ggplot(fp_data,
               aes(x = loci_set, y = population, fill = prop_gt_1)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.4f\n(n=%d)", prop_gt_1, n_gt_1)),
              size = 3) +
    facet_wrap(~ tested_relationship) +
    scale_fill_gradient(low = "white", high = "red",
                        name = "Prop > 1",
                        limits = c(0, NA)) +
    labs(title = "False Positive Rate: Unrelated Pairs (Correct Population)",
         subtitle = "Proportion with combined LR > 1",
         x = "Loci Set", y = "True Population") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold"))

  ggsave(file.path(input_dir, "fig1_unrelated_fp_heatmap.png"),
         p1, width = 12, height = 6, dpi = 150)
}

# ==========================================================================
# Figure 2: True Positive rates across loci sets (correct pop)
# ==========================================================================
cat("Generating Figure 2: True positive rates...\n")

tp_data <- dt_pop[classification == "True Positive" & is_correct_pop == TRUE]

if (nrow(tp_data) > 0) {
  tp_long <- tp_data %>%
    select(population, known_relationship, tested_relationship, loci_set,
           prop_gt_1, prop_gt_10, prop_gt_100, prop_gt_1000) %>%
    pivot_longer(cols = starts_with("prop_gt_"),
                 names_to = "threshold",
                 values_to = "proportion") %>%
    mutate(threshold = factor(threshold,
                              levels = c("prop_gt_1", "prop_gt_10",
                                         "prop_gt_100", "prop_gt_1000"),
                              labels = c("LR > 1", "LR > 10",
                                         "LR > 100", "LR > 1000")))

  p2 <- ggplot(tp_long,
               aes(x = loci_set, y = proportion,
                   color = threshold, group = threshold)) +
    geom_point(size = 2) +
    geom_line() +
    facet_grid(population ~ tested_relationship) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    scale_color_brewer(palette = "Set1", name = "Threshold") +
    labs(title = "True Positive Rate by Loci Set (Correct Population)",
         subtitle = "Proportion of truly related pairs exceeding LR thresholds",
         x = "Loci Set", y = "Proportion Exceeding Threshold") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold"))

  ggsave(file.path(input_dir, "fig2_true_positive_rates.png"),
         p2, width = 12, height = 10, dpi = 150)
}

# ==========================================================================
# Figure 3: Cross-relationship FP rates (correct pop)
# ==========================================================================
cat("Generating Figure 3: Cross-relationship FP...\n")

cross_data <- dt_pop[classification == "Related FP" & is_correct_pop == TRUE]

if (nrow(cross_data) > 0) {
  p3 <- ggplot(cross_data,
               aes(x = loci_set, y = prop_gt_100,
                   color = known_relationship, group = known_relationship)) +
    geom_point(size = 2) +
    geom_line() +
    facet_grid(population ~ tested_relationship) +
    scale_y_continuous(labels = scales::percent) +
    scale_color_brewer(palette = "Set2", name = "True Relationship") +
    labs(title = "Cross-Relationship False Positive Rate (Correct Population)",
         subtitle = "Proportion of misclassified related pairs with LR > 100",
         x = "Loci Set", y = "Proportion with LR > 100") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold"))

  ggsave(file.path(input_dir, "fig3_cross_relationship_fp.png"),
         p3, width = 12, height = 10, dpi = 150)
}

# ==========================================================================
# Figure 4: Correct vs Wrong Population comparison
# ==========================================================================
cat("Generating Figure 4: Population mismatch effect...\n")

# Compare correct pop vs wrong pop for TP cases
pop_cmp <- dt_pop[classification == "True Positive"]

if (nrow(pop_cmp) > 0 && any(pop_cmp$is_correct_pop == FALSE)) {
  p4 <- ggplot(pop_cmp,
               aes(x = loci_set, y = mean_log10_LR,
                   color = is_correct_pop, group = interaction(tested_population, is_correct_pop))) +
    geom_point(size = 2, alpha = 0.7) +
    geom_line(alpha = 0.5) +
    facet_grid(population ~ tested_relationship) +
    scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "red"),
                       labels = c("TRUE" = "Correct Pop", "FALSE" = "Wrong Pop"),
                       name = "Population Match") +
    labs(title = "Effect of Population Mismatch on LR (True Positive Pairs)",
         subtitle = "Mean log10(LR) for correct vs wrong population frequency tables",
         x = "Loci Set", y = "Mean log10(LR)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold"))

  ggsave(file.path(input_dir, "fig4_population_mismatch.png"),
         p4, width = 12, height = 10, dpi = 150)
} else {
  cat("  Skipped (no wrong-population data available)\n")
}

# ==========================================================================
# Figure 5: Overall classification bar chart
# ==========================================================================
cat("Generating Figure 5: Classification overview...\n")

overview <- dt_pop[is_correct_pop == TRUE & loci_set == "Core 13"]

if (nrow(overview) > 0) {
  overview_long <- overview %>%
    select(population, known_relationship, tested_relationship,
           classification, prop_gt_1, prop_gt_100, prop_gt_1000) %>%
    pivot_longer(cols = starts_with("prop_gt_"),
                 names_to = "threshold",
                 values_to = "proportion") %>%
    mutate(threshold = factor(threshold,
                              levels = c("prop_gt_1", "prop_gt_100", "prop_gt_1000"),
                              labels = c("LR > 1", "LR > 100", "LR > 1000")))

  p5 <- ggplot(overview_long,
               aes(x = known_relationship, y = proportion, fill = classification)) +
    geom_col(position = "dodge") +
    facet_grid(threshold ~ tested_relationship, scales = "free_y") +
    scale_fill_manual(values = c("True Positive" = "#2ca02c",
                                 "Related FP" = "#ff7f0e",
                                 "Unrelated FP" = "#d62728"),
                      name = "Classification") +
    scale_y_continuous(labels = scales::percent) +
    labs(title = "Classification Overview: Core 13 Loci (Correct Population)",
         subtitle = "Proportion exceeding LR threshold by true relationship",
         x = "True Relationship", y = "Proportion") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold"))

  ggsave(file.path(input_dir, "fig5_classification_overview.png"),
         p5, width = 14, height = 10, dpi = 150)
}

# ==========================================================================
# Figure 6: LR=0 diagnostic — where are zeros happening?
# ==========================================================================
cat("Generating Figure 6: Zero LR diagnostic...\n")

zero_data <- dt_pop[is_correct_pop == TRUE,
                    .(population, known_relationship, tested_relationship,
                      loci_set, n_zero, n_pairs,
                      prop_zero = ifelse(n_pairs > 0, n_zero / n_pairs, 0))]

if (any(zero_data$n_zero > 0)) {
  p6 <- ggplot(zero_data[n_zero > 0],
               aes(x = loci_set, y = known_relationship,
                   fill = prop_zero)) +
    geom_tile(color = "white") +
    geom_text(aes(label = n_zero), size = 3) +
    facet_grid(population ~ tested_relationship) +
    scale_fill_gradient(low = "lightyellow", high = "darkred",
                        name = "Prop Zero") +
    labs(title = "Diagnostic: Combined LR = 0 (Correct Population)",
         subtitle = "Count of pairs with LR exactly zero (at least one locus with 0 IBS sharing)",
         x = "Loci Set", y = "True Relationship") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold"))

  ggsave(file.path(input_dir, "fig6_zero_lr_diagnostic.png"),
         p6, width = 12, height = 10, dpi = 150)
} else {
  cat("  No zero LR values found, skipping.\n")
}

cat("\nAll figures saved to:", input_dir, "\n")
cat("=== Done ===\n")
