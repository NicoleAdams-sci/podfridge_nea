#!/usr/bin/env Rscript

# =============================================================================
# Prepare Intermediate Files for Publication Plotting Scripts
# =============================================================================
# Loads combined_LR_all.rds once and writes three small aggregated CSVs
# consumed by the downstream publication plotting scripts. This avoids
# reloading the full 5M-row dataset on every figure iteration.
#
# Outputs (written to <output_dir>/):
#   proportions_with_classification.csv  → plots_cutoffs_publication.R
#   mismatched_pop_robustness.csv        → plots_mismatched_population.R (main figure + table)
#   mismatched_pop_heatmap.csv           → plots_mismatched_population.R (supplement heatmap)
#
# Usage:
#   Rscript code/prepare_combined_lr_intermediates.R <input_dir> [output_dir]
#
#   input_dir   Full path to directory containing combined_LR_all.rds
#               (e.g., output/lr_analysis_20260130)
#   output_dir  Where to write intermediate CSVs
#               (default: <input_dir>, alongside other analysis outputs)
# =============================================================================

suppressMessages(suppressWarnings({
  library(tidyverse)
  library(data.table)
}))

source("code/module9_combinedLR_stats_functions.R")

log_message <- function(msg) cat(paste0("[", Sys.time(), "] ", msg, "\n"))

# --- Argument parsing ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript code/prepare_combined_lr_intermediates.R <input_dir> [output_dir]")
}

input_dir  <- args[1]
output_dir <- if (length(args) >= 2) args[2] else input_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

log_message(paste("Input directory: ", input_dir))
log_message(paste("Output directory:", output_dir))


# =============================================================================
# SECTION 1: SHARED CONSTANTS
# =============================================================================

relationship_order <- c("parent_child", "full_siblings", "half_siblings",
                        "cousins", "second_cousins", "unrelated")

loci_set_order <- c("core_13", "identifiler_15", "expanded_20",
                    "supplementary", "autosomal_29")

population_order <- c("AfAm", "Asian", "Cauc", "Hispanic", "all")


# =============================================================================
# SECTION 2: LOAD RAW DATA
# =============================================================================

log_message("Loading combined_LR_all.rds...")

all_combined_file <- file.path(input_dir, "combined_LR_all.rds")
if (!file.exists(all_combined_file)) {
  stop(sprintf("Data file not found: %s", all_combined_file))
}

all_combined <- readRDS(all_combined_file)
log_message(sprintf("Loaded %s rows.", format(nrow(all_combined), big.mark = ",")))


# =============================================================================
# SECTION 3: DIAGNOSTICS (pre-filter)
# Runs once here so quality checks don't need to be repeated in each plotting
# script. Results are written to the SLURM log for post-hoc inspection.
# =============================================================================

log_message("Running pre-filter diagnostics...")

all_combined <- all_combined %>% mutate(combined_LR = as.numeric(combined_LR))

n_total_raw  <- nrow(all_combined)
n_zero_raw   <- sum(all_combined$combined_LR == 0, na.rm = TRUE)
n_na_raw     <- sum(is.na(all_combined$combined_LR))
n_neginf_raw <- sum(!is.finite(log10(suppressWarnings(
                  all_combined$combined_LR))), na.rm = TRUE)

cat(sprintf("\n--- LR Diagnostic (pre-filter, n = %s) ---\n",
            format(n_total_raw, big.mark = ",")))
cat(sprintf("  Zero LRs:              %d (%.4f%%)\n",
            n_zero_raw,   100 * n_zero_raw   / n_total_raw))
cat(sprintf("  NA LRs:                %d (%.4f%%)\n",
            n_na_raw,     100 * n_na_raw     / n_total_raw))
cat(sprintf("  Non-finite log10(LR):  %d (%.4f%%)\n\n",
            n_neginf_raw, 100 * n_neginf_raw / n_total_raw))

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

cat("--- Exclusion rate by population match (pre-filter) ---\n")
print(as.data.frame(excl_by_match))

zero_by_relationship <- all_combined %>%
  filter(combined_LR == 0) %>%
  count(known_relationship, tested_relationship, name = "n_zeros") %>%
  arrange(desc(n_zeros)) %>%
  mutate(pct_of_all_zeros = 100 * n_zeros / n_zero_raw)

cat("\n--- Zero LR breakdown by known vs. tested relationship (pre-filter) ---\n")
print(as.data.frame(zero_by_relationship))
cat("\n")


# =============================================================================
# SECTION 4: PREPROCESS
# =============================================================================

log_message("Preprocessing...")

all_combined <- all_combined %>%
  mutate(
    known_relationship  = factor(known_relationship,  levels = relationship_order),
    tested_relationship = factor(tested_relationship, levels = relationship_order),
    loci_set            = factor(loci_set,            levels = loci_set_order),
    population          = factor(population,          levels = population_order),
    pop_match           = tested_population == population,
    pop_match_label     = if_else(pop_match, "Matched", "Mismatched"),
    rel_match           = as.character(known_relationship) == as.character(tested_relationship),
    log10_LR            = suppressWarnings(log10(combined_LR))
  ) %>%
  filter(!is.na(combined_LR))

n_total_clean <- nrow(all_combined)
n_zero_clean  <- sum(all_combined$combined_LR == 0, na.rm = TRUE)
cat(sprintf("\n--- LR Diagnostic (post-filter, n = %s) ---\n",
            format(n_total_clean, big.mark = ",")))
cat(sprintf("  Zero LRs retained: %d (%.4f%%)\n\n",
            n_zero_clean, 100 * n_zero_clean / n_total_clean))

log_message(sprintf("After preprocessing: %s rows.", format(nrow(all_combined), big.mark = ",")))


# =============================================================================
# SECTION 5: proportions_with_classification.csv
# Consumed by: plots_cutoffs_publication.R
#
# Calculates empirical cutoffs from unrelated pairs, then the proportion of
# each relationship x loci x population group exceeding each threshold.
# Adds classification labels and convenience prop_LR_gt_* column aliases.
# =============================================================================

log_message("Building proportions_with_classification...")

cutoffs <- calculate_cutoffs(all_combined, c(1, 0.1, 0.01))
proportions_exceeding <- calculate_proportions_exceeding_cutoffs(all_combined, cutoffs)

proportions_with_classification <- proportions_exceeding %>%
  mutate(
    classification = case_when(
      as.character(known_relationship) == as.character(tested_relationship) ~ "True Positive",
      known_relationship == "unrelated"                                     ~ "Unrelated FP",
      TRUE                                                                  ~ "Related FP"
    ),
    n_loci = case_when(
      loci_set == "core_13"        ~ 13,
      loci_set == "identifiler_15" ~ 15,
      loci_set == "expanded_20"    ~ 20,
      loci_set == "supplementary"  ~ 23,
      loci_set == "autosomal_29"   ~ 29
    ),
    prop_LR_gt_1    = proportion_exceeding_fixed,
    prop_LR_gt_10   = proportion_exceeding_1,
    prop_LR_gt_100  = proportion_exceeding_0_1,
    prop_LR_gt_1000 = proportion_exceeding_0_01
  )

out1 <- file.path(output_dir, "proportions_with_classification.csv")
fwrite(proportions_with_classification, out1)
log_message(sprintf("Wrote %d rows -> %s", nrow(proportions_with_classification), out1))


# =============================================================================
# SECTION 6: mismatched_pop_robustness.csv
# Consumed by: plots_mismatched_population.R (main figure + summary table)
#
# True positive pairs only (relationship correctly tested). Summarizes median
# and mean log10(LR) per population x relationship x loci x match status cell.
# =============================================================================

log_message("Building mismatched_pop_robustness...")

mismatched_pop_robustness <- all_combined %>%
  filter(
    rel_match == TRUE,
    known_relationship %in% c("parent_child", "full_siblings")
  ) %>%
  group_by(population, known_relationship, loci_set, pop_match_label, tested_population) %>%
  summarize(
    median_log10_LR = median(log10_LR, na.rm = TRUE),
    mean_log10_LR   = mean(log10_LR[is.finite(log10_LR)], na.rm = TRUE),
    se_log10_LR     = sd(log10_LR[is.finite(log10_LR)], na.rm = TRUE) /
                        sqrt(sum(combined_LR > 0)),
    n_pairs         = n(),
    n_zero          = sum(combined_LR == 0),
    .groups = "drop"
  ) %>%
  mutate(
    n_loci = case_when(
      loci_set == "core_13"        ~ 13,
      loci_set == "identifiler_15" ~ 15,
      loci_set == "expanded_20"    ~ 20,
      loci_set == "supplementary"  ~ 23,
      loci_set == "autosomal_29"   ~ 29
    )
  )

out2 <- file.path(output_dir, "mismatched_pop_robustness.csv")
fwrite(mismatched_pop_robustness, out2)
log_message(sprintf("Wrote %d rows -> %s", nrow(mismatched_pop_robustness), out2))


# =============================================================================
# SECTION 7: mismatched_pop_heatmap.csv
# Consumed by: plots_mismatched_population.R (supplement heatmap)
#
# For each true population x tested population x relationship x loci cell,
# computes the median log10(LR_wrong / LR_correct) across all pairs.
# Off-diagonal rows are joined by pair ID; diagonal rows (pop_match == TRUE)
# have log10_ratio = 0 by construction and are added directly.
# =============================================================================

log_message("Building mismatched_pop_heatmap...")

# Off-diagonal: pairs where wrong population frequencies were used
correct_lrs <- all_combined %>%
  filter(
    pop_match == TRUE, rel_match == TRUE,
    known_relationship %in% c("parent_child", "full_siblings")
  ) %>%
  select(batch_id, pair_id, population, known_relationship,
         loci_set, tested_relationship,
         correct_log10_LR = log10_LR)

wrong_lrs <- all_combined %>%
  filter(
    pop_match == FALSE, rel_match == TRUE,
    known_relationship %in% c("parent_child", "full_siblings")
  ) %>%
  select(batch_id, pair_id, population, known_relationship,
         loci_set, tested_relationship, tested_population,
         wrong_log10_LR = log10_LR)

off_diagonal <- wrong_lrs %>%
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
  mutate(is_diagonal = FALSE)

# Diagonal: matched frequencies — LR_wrong / LR_correct = 1, so log10_ratio = 0
on_diagonal <- all_combined %>%
  filter(
    pop_match == TRUE, rel_match == TRUE,
    known_relationship %in% c("parent_child", "full_siblings")
  ) %>%
  group_by(population, tested_population, known_relationship, loci_set) %>%
  summarize(
    median_log10_ratio = 0,
    n_pairs            = n(),
    .groups = "drop"
  ) %>%
  mutate(is_diagonal = TRUE)

mismatched_pop_heatmap <- bind_rows(off_diagonal, on_diagonal)

out3 <- file.path(output_dir, "mismatched_pop_heatmap.csv")
fwrite(mismatched_pop_heatmap, out3)
log_message(sprintf("Wrote %d rows -> %s", nrow(mismatched_pop_heatmap), out3))


# =============================================================================
# SUMMARY
# =============================================================================

cat("\nINTERMEDIATE FILES WRITTEN:\n")
cat(sprintf("  %s\n", out1))
cat(sprintf("  %s\n", out2))
cat(sprintf("  %s\n", out3))
cat(sprintf("\nAll files saved to: %s\n", output_dir))
log_message("Done.")
