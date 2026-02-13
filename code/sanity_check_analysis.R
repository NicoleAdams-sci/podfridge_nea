#!/usr/bin/env Rscript
# ==============================================================================
# sanity_check_analysis.R
# Comprehensive sanity check: false positive/negative rates across ALL combos
#
# Produces a single CSV spreadsheet with counts and proportions for:
#   - Tested hypotheses: parent_child, full_siblings
#   - Known relationships: all 6
#   - Correct vs wrong population assumptions
#   - All 5 loci sets
#   - LR thresholds: 1, 10, 100, 1000
#
# Usage: Rscript code/sanity_check_analysis.R [OUTPUT_DIR] [COMBINED_LR_DIR]
# ==============================================================================

library(data.table)
library(dplyr)

cat("=== Sanity Check Analysis ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# --------------------------------------------------------------------------
# 0. Parse arguments
# --------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

output_dir <- if (length(args) >= 1) args[1] else "output/sanity_check"
combined_lr_dir <- if (length(args) >= 2) args[2] else "output/combined_LR"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cat("Output dir:", output_dir, "\n")
cat("Reading combined_LR files from:", combined_lr_dir, "\n\n")

# --------------------------------------------------------------------------
# 1. Read all combined_LR files
# --------------------------------------------------------------------------
all_files <- list.files(combined_lr_dir,
                        pattern = "^combined_LR_.*\\.csv$",
                        full.names = TRUE)

if (length(all_files) == 0) {
  stop("No combined_LR CSV files found in ", combined_lr_dir)
}

cat("Found", length(all_files), "combined_LR files\n")

# Read efficiently with data.table
dt_list <- lapply(all_files, function(f) {
  tryCatch(fread(f), error = function(e) {
    cat("  Warning: Failed to read", f, "-", e$message, "\n")
    NULL
  })
})
dt_list <- dt_list[!sapply(dt_list, is.null)]
all_combined <- rbindlist(dt_list, fill = TRUE)
rm(dt_list)
gc()

cat("Total rows loaded:", nrow(all_combined), "\n")
all_combined[, combined_LR := as.numeric(combined_LR)]

# --------------------------------------------------------------------------
# 2. Diagnostic: check what we have
# --------------------------------------------------------------------------
cat("\n--- Data Summary ---\n")
cat("Known relationships:", paste(sort(unique(all_combined$known_relationship)), collapse = ", "), "\n")
cat("Tested relationships:", paste(sort(unique(all_combined$tested_relationship)), collapse = ", "), "\n")
cat("Populations (true):", paste(sort(unique(all_combined$population)), collapse = ", "), "\n")
cat("Tested populations:", paste(sort(unique(all_combined$tested_population)), collapse = ", "), "\n")
cat("Loci sets:", paste(sort(unique(all_combined$loci_set)), collapse = ", "), "\n")

# Count unique pairs per known_relationship x population
pair_counts <- all_combined[, .(n_unique_pairs = uniqueN(paste(batch_id, pair_id))),
                            by = .(population, known_relationship)]
cat("\n--- Unique pairs per population x relationship ---\n")
print(pair_counts[order(population, known_relationship)])

# --------------------------------------------------------------------------
# 3. Compute FALSE POSITIVE / FALSE NEGATIVE table
#    Focus on tested_relationship = parent_child and full_siblings
# --------------------------------------------------------------------------
tested_hyps <- c("parent_child", "full_siblings")

# Filter to tested hypotheses of interest
dt_focus <- all_combined[tested_relationship %in% tested_hyps]

cat("\nFiltered to tested hypotheses:", paste(tested_hyps, collapse = ", "), "\n")
cat("Rows after filter:", nrow(dt_focus), "\n")

# Compute counts and proportions at standard LR thresholds
thresholds <- c(1, 10, 100, 1000)

results <- dt_focus[, {
  n_total <- .N
  # Number of unique pairs (sanity check)
  n_pairs <- uniqueN(paste(batch_id, pair_id))

  # Count how many exceed each threshold
  n_gt_1    <- sum(combined_LR > 1, na.rm = TRUE)
  n_gt_10   <- sum(combined_LR > 10, na.rm = TRUE)
  n_gt_100  <- sum(combined_LR > 100, na.rm = TRUE)
  n_gt_1000 <- sum(combined_LR > 1000, na.rm = TRUE)

  # Count NAs and zeros
  n_na    <- sum(is.na(combined_LR))
  n_zero  <- sum(combined_LR == 0, na.rm = TRUE)
  n_inf   <- sum(is.infinite(combined_LR))

  # Proportions
  prop_gt_1    <- n_gt_1 / n_total
  prop_gt_10   <- n_gt_10 / n_total
  prop_gt_100  <- n_gt_100 / n_total
  prop_gt_1000 <- n_gt_1000 / n_total

  # Summary LR stats
  mean_log10_LR <- mean(log10(pmax(combined_LR, 1e-300)), na.rm = TRUE)
  median_log10_LR <- median(log10(pmax(combined_LR, 1e-300)), na.rm = TRUE)

  list(
    n_rows     = n_total,
    n_pairs    = n_pairs,
    n_na       = n_na,
    n_zero     = n_zero,
    n_inf      = n_inf,
    n_gt_1     = n_gt_1,
    n_gt_10    = n_gt_10,
    n_gt_100   = n_gt_100,
    n_gt_1000  = n_gt_1000,
    prop_gt_1  = round(prop_gt_1, 6),
    prop_gt_10 = round(prop_gt_10, 6),
    prop_gt_100  = round(prop_gt_100, 6),
    prop_gt_1000 = round(prop_gt_1000, 6),
    mean_log10_LR   = round(mean_log10_LR, 3),
    median_log10_LR = round(median_log10_LR, 3)
  )
}, by = .(population, known_relationship, tested_relationship,
          loci_set, tested_population, is_correct_pop)]

# --------------------------------------------------------------------------
# 4. Add classification labels
# --------------------------------------------------------------------------
results[, classification := fcase(
  known_relationship == tested_relationship, "True Positive",
  known_relationship == "unrelated",         "Unrelated FP",
  default = "Related FP"
)]

# For false negatives: when known is related but LR < threshold
# (We can read these from prop columns: FN rate = 1 - prop_gt_threshold for TP rows)

# Sort for readability
setorder(results, population, tested_relationship, loci_set,
         is_correct_pop, tested_population, known_relationship)

# --------------------------------------------------------------------------
# 5. Save comprehensive spreadsheet
# --------------------------------------------------------------------------
out_file <- file.path(output_dir, "sanity_check_all_combinations.csv")
fwrite(results, out_file)
cat("\nWrote comprehensive results to:", out_file, "\n")
cat("Total rows in output:", nrow(results), "\n")

# --------------------------------------------------------------------------
# 6. Create focused summary tables
# --------------------------------------------------------------------------

# 6a. FALSE POSITIVE SUMMARY (unrelated pairs tested as PC or FS)
#     This is the key sanity check: how many unrelated pass?
fp_summary <- results[known_relationship == "unrelated",
                      .(population, tested_relationship, loci_set,
                        tested_population, is_correct_pop,
                        n_pairs, n_gt_1, n_gt_10, n_gt_100, n_gt_1000,
                        prop_gt_1, prop_gt_10, prop_gt_100, prop_gt_1000,
                        mean_log10_LR)]

setorder(fp_summary, tested_relationship, loci_set, population, is_correct_pop)
fp_file <- file.path(output_dir, "sanity_check_false_positives_unrelated.csv")
fwrite(fp_summary, fp_file)
cat("Wrote unrelated FP summary to:", fp_file, "\n")

# 6b. FALSE NEGATIVE SUMMARY (related pairs that DON'T pass)
#     Only for true positive rows (known == tested)
fn_summary <- results[known_relationship == tested_relationship,
                      .(population, known_relationship, tested_relationship,
                        loci_set, tested_population, is_correct_pop,
                        n_pairs,
                        fn_rate_t1    = round(1 - prop_gt_1, 6),
                        fn_rate_t10   = round(1 - prop_gt_10, 6),
                        fn_rate_t100  = round(1 - prop_gt_100, 6),
                        fn_rate_t1000 = round(1 - prop_gt_1000, 6),
                        mean_log10_LR)]

setorder(fn_summary, known_relationship, loci_set, population, is_correct_pop)
fn_file <- file.path(output_dir, "sanity_check_false_negatives.csv")
fwrite(fn_summary, fn_file)
cat("Wrote FN summary to:", fn_file, "\n")

# 6c. CROSS-RELATIONSHIP FP: related pairs tested under wrong hypothesis
cross_fp <- results[known_relationship != tested_relationship &
                    known_relationship != "unrelated",
                    .(population, known_relationship, tested_relationship,
                      loci_set, tested_population, is_correct_pop,
                      n_pairs, n_gt_1, n_gt_10, n_gt_100, n_gt_1000,
                      prop_gt_1, prop_gt_10, prop_gt_100, prop_gt_1000,
                      classification)]

setorder(cross_fp, tested_relationship, known_relationship, loci_set,
         population, is_correct_pop)
cross_fp_file <- file.path(output_dir, "sanity_check_cross_relationship_fp.csv")
fwrite(cross_fp, cross_fp_file)
cat("Wrote cross-relationship FP to:", cross_fp_file, "\n")

# 6d. POPULATION MISMATCH SUMMARY: correct pop vs wrong pop comparison
pop_compare <- results[, .(population, known_relationship, tested_relationship,
                           loci_set, tested_population, is_correct_pop,
                           n_pairs, prop_gt_1, prop_gt_100, prop_gt_1000,
                           mean_log10_LR)]

# Pivot to wide: correct vs wrong side by side
pop_correct <- pop_compare[is_correct_pop == TRUE,
                           .(population, known_relationship, tested_relationship,
                             loci_set, tested_population,
                             n_pairs_correct = n_pairs,
                             prop_gt_1_correct = prop_gt_1,
                             prop_gt_100_correct = prop_gt_100,
                             prop_gt_1000_correct = prop_gt_1000,
                             mean_log10_LR_correct = mean_log10_LR)]

pop_wrong <- pop_compare[is_correct_pop == FALSE,
                         .(population, known_relationship, tested_relationship,
                           loci_set, tested_population,
                           n_pairs_wrong = n_pairs,
                           prop_gt_1_wrong = prop_gt_1,
                           prop_gt_100_wrong = prop_gt_100,
                           prop_gt_1000_wrong = prop_gt_1000,
                           mean_log10_LR_wrong = mean_log10_LR)]

# Save both
pop_compare_file <- file.path(output_dir, "sanity_check_population_comparison.csv")
fwrite(pop_compare, pop_compare_file)
cat("Wrote population comparison to:", pop_compare_file, "\n")

# --------------------------------------------------------------------------
# 7. Print key sanity check results to console
# --------------------------------------------------------------------------
cat("\n")
cat("===============================================================\n")
cat("                     KEY SANITY CHECKS                         \n")
cat("===============================================================\n\n")

# Check 1: Unrelated pairs tested as parent_child with correct population
cat("--- Unrelated pairs tested as PARENT_CHILD (correct pop) ---\n")
check1 <- results[known_relationship == "unrelated" &
                  tested_relationship == "parent_child" &
                  is_correct_pop == TRUE]
if (nrow(check1) > 0) {
  for (i in 1:nrow(check1)) {
    r <- check1[i]
    cat(sprintf("  Pop=%s, Loci=%s, Tested_pop=%s: n=%d, FP_rate(>1)=%.6f, FP_rate(>100)=%.6f, FP_rate(>1000)=%.6f\n",
                r$population, r$loci_set, r$tested_population,
                r$n_pairs, r$prop_gt_1, r$prop_gt_100, r$prop_gt_1000))
  }
} else {
  cat("  No data found!\n")
}

cat("\n--- Unrelated pairs tested as FULL_SIBLINGS (correct pop) ---\n")
check2 <- results[known_relationship == "unrelated" &
                  tested_relationship == "full_siblings" &
                  is_correct_pop == TRUE]
if (nrow(check2) > 0) {
  for (i in 1:nrow(check2)) {
    r <- check2[i]
    cat(sprintf("  Pop=%s, Loci=%s, Tested_pop=%s: n=%d, FP_rate(>1)=%.6f, FP_rate(>100)=%.6f, FP_rate(>1000)=%.6f\n",
                r$population, r$loci_set, r$tested_population,
                r$n_pairs, r$prop_gt_1, r$prop_gt_100, r$prop_gt_1000))
  }
} else {
  cat("  No data found!\n")
}

cat("\n--- Parent-child pairs tested as PARENT_CHILD (correct pop, FN check) ---\n")
check3 <- results[known_relationship == "parent_child" &
                  tested_relationship == "parent_child" &
                  is_correct_pop == TRUE]
if (nrow(check3) > 0) {
  for (i in 1:nrow(check3)) {
    r <- check3[i]
    cat(sprintf("  Pop=%s, Loci=%s: n=%d, TP_rate(>1)=%.4f, TP_rate(>100)=%.4f, TP_rate(>1000)=%.4f\n",
                r$population, r$loci_set,
                r$n_pairs, r$prop_gt_1, r$prop_gt_100, r$prop_gt_1000))
  }
} else {
  cat("  No data found!\n")
}

cat("\n--- Full-sib pairs tested as FULL_SIBLINGS (correct pop, FN check) ---\n")
check4 <- results[known_relationship == "full_siblings" &
                  tested_relationship == "full_siblings" &
                  is_correct_pop == TRUE]
if (nrow(check4) > 0) {
  for (i in 1:nrow(check4)) {
    r <- check4[i]
    cat(sprintf("  Pop=%s, Loci=%s: n=%d, TP_rate(>1)=%.4f, TP_rate(>100)=%.4f, TP_rate(>1000)=%.4f\n",
                r$population, r$loci_set,
                r$n_pairs, r$prop_gt_1, r$prop_gt_100, r$prop_gt_1000))
  }
} else {
  cat("  No data found!\n")
}

# Check for any zeros in combined_LR for parent_child tested as parent_child
cat("\n--- DIAGNOSTIC: LR=0 counts (indicates allele sharing issues) ---\n")
zero_check <- all_combined[combined_LR == 0 & !is.na(combined_LR),
                           .N, by = .(known_relationship, tested_relationship, loci_set)]
if (nrow(zero_check) > 0) {
  setorder(zero_check, known_relationship, tested_relationship, loci_set)
  print(zero_check)
} else {
  cat("  No LR=0 values found (good!)\n")
}

cat("\n--- DIAGNOSTIC: LR=NA counts ---\n")
na_check <- all_combined[is.na(combined_LR),
                         .N, by = .(known_relationship, tested_relationship, loci_set)]
if (nrow(na_check) > 0) {
  setorder(na_check, known_relationship, tested_relationship, loci_set)
  print(na_check)
} else {
  cat("  No LR=NA values found\n")
}

cat("\nEnd time:", format(Sys.time()), "\n")
cat("=== Done ===\n")
