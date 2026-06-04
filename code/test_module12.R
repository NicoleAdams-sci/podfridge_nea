##### Module 12 Standalone Test #####
# USAGE: Rscript code/test_module12.R <LR_COMBINED_FILE>
#
# Example:
#   Rscript code/test_module12.R output/focal_ranking_test/ranking_lrs_all_replicates_second_cousins_20260604_172400.csv
#
# Purpose:
#   Tests module12_ranking_outcome_recorder.R by:
#     1. Loading the combined LR results from the specified file
#     2. Running record_ranking_outcomes() on the combined data
#     3. Checking output structure and sanity of results
#
# Prerequisites:
#   Run focal_combine_ranking_lrs.R first to generate the combined file.

# ------------------------------------------------------------------------------
# 0. Parse arguments
# ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript test_module12.R <LR_COMBINED_FILE>")
}

LR_COMBINED_FILE <- args[1]

if (!file.exists(LR_COMBINED_FILE))
  stop(paste("File not found:", LR_COMBINED_FILE))

cat("=== Module 12 Standalone Test ===\n")
cat("Started at:", format(Sys.time()), "\n\n")

# ------------------------------------------------------------------------------
# 0. Load packages and source module 12
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

source("code/module12_ranking_outcome_recorder.R")

TOP_N      <- 200
OUTPUT_DIR <- "output/focal_ranking_test"

# ------------------------------------------------------------------------------
# STEP 1: Load combined LR results
# ------------------------------------------------------------------------------

cat("--- STEP 1: Load combined LR results ---\n")

cat(sprintf("Loading: %s\n", basename(LR_COMBINED_FILE)))

lr_results <- fread(LR_COMBINED_FILE)

# Ensure combined_LR is numeric — fread may read as character
# if scientific notation values are mixed with other types in the CSV
lr_results[, combined_LR := as.numeric(combined_LR)]

cat(sprintf("Rows loaded              : %s\n", format(nrow(lr_results), big.mark = ",")))
cat(sprintf("Unique focal_ids         : %s\n", paste(sort(unique(lr_results$focal_id)), collapse = ", ")))
cat(sprintf("Loci sets                : %s\n", paste(unique(lr_results$loci_set), collapse = ", ")))
cat(sprintf("Tested relationships     : %s\n", paste(unique(lr_results$tested_relationship), collapse = ", ")))
cat(sprintf("Tested populations       : %s\n", paste(unique(lr_results$tested_population), collapse = ", ")))
cat(sprintf("is_true_relative present : %s\n\n", ifelse("is_true_relative" %in% names(lr_results), "YES", "NO")))

# Quick check — confirm exactly one true relative per focal x loci_set x tested_rel
true_rel_counts <- lr_results |>
  filter(is_true_relative == TRUE) |>
  group_by(focal_id, loci_set, tested_relationship) |>
  summarize(n = n(), .groups = "drop") |>
  pull(n) |> unique()

cat(sprintf("True relatives per focal x loci_set x tested_rel: %s (expected: all 1)\n\n",
            paste(true_rel_counts, collapse = ", ")))

# ------------------------------------------------------------------------------
# STEP 2: Run record_ranking_outcomes() (module 12)
# ------------------------------------------------------------------------------

cat("--- STEP 2: Run record_ranking_outcomes() ---\n\n")

results <- record_ranking_outcomes(
  lr_results   = lr_results,
  top_n        = TOP_N,
  output_dir   = OUTPUT_DIR,
  save_results = TRUE
)

outcomes <- results$outcomes
summary  <- results$summary

# ------------------------------------------------------------------------------
# STEP 3: Checks
# ------------------------------------------------------------------------------

cat("\n--- STEP 3: Output checks ---\n")

# Expected: one row per focal_id x loci_set x tested_relationship x tested_population
n_focal     <- length(unique(lr_results$focal_id))
n_loci_sets <- length(unique(lr_results$loci_set))
n_tested    <- length(unique(lr_results$tested_relationship))
n_pops      <- length(unique(lr_results$tested_population))
expected_outcome_rows <- n_focal * n_loci_sets * n_tested * n_pops

cat(sprintf("Outcomes rows   : %d (expected: %d)\n", nrow(outcomes), expected_outcome_rows))
cat(sprintf("Summary rows    : %d\n", nrow(summary)))
cat(sprintf("Rank range      : %g - %g\n", min(outcomes$rank), max(outcomes$rank)))
cat(sprintf("NA ranks        : %d (expected: 0)\n", sum(is.na(outcomes$rank))))
cat(sprintf("NA in_top_n     : %d (expected: 0)\n", sum(is.na(outcomes$in_top_n))))

# ------------------------------------------------------------------------------
# STEP 4: Print full summary table
# ------------------------------------------------------------------------------

cat("\n--- STEP 4: Full summary table ---\n")
print(as.data.frame(summary |>
  select(loci_set, tested_relationship, known_relationship,
         n_replicates, prop_in_top_n, median_rank, min_rank, max_rank) |>
  arrange(loci_set, tested_relationship)))

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------

cat("\n=== Module 12 Test Summary ===\n")

checks <- c(
  "Correct outcomes row count" = nrow(outcomes) == expected_outcome_rows,
  "No NA ranks"                = sum(is.na(outcomes$rank)) == 0,
  "No NA in_top_n"             = sum(is.na(outcomes$in_top_n)) == 0,
  "Ranks within database size" = max(outcomes$rank) <= 1001,
  "Summary rows match groups"  = nrow(summary) == n_loci_sets * n_tested * n_pops
)

for (check_name in names(checks)) {
  status <- ifelse(isTRUE(checks[[check_name]]), "PASS", "FAIL")
  cat(sprintf("  [%s] %s\n", status, check_name))
}

if (all(sapply(checks, isTRUE))) {
  cat("\nSUCCESS: All checks passed. Module 12 is working correctly.\n")
} else {
  cat("\nWARNING: Some checks failed. Review output above.\n")
}

cat("Completed at:", format(Sys.time()), "\n")
