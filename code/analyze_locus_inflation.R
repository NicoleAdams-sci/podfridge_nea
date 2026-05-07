#!/usr/bin/env Rscript
# =============================================================================
# analyze_locus_inflation.R
#
# Purpose: Identify which loci drive LR inflation when using mismatched
#          population allele frequencies.
#
# Two outputs:
#   1. locus_inflation_summary.csv
#      Per-locus median log10(LR_wrong / LR_correct) for each true population
#      x tested population combination. Restricted to true positives
#      (known_relationship == tested_relationship). Ranks loci by their
#      contribution to inflation, with focus on AfAm <-> Asian mismatch.
#
#   2. locus_heterozygosity_summary.csv
#      Per-locus observed heterozygosity per true population, computed from
#      simulated genotypes (focal_A1, focal_A2). Used to assess whether loci
#      driving the most inflation are also those with the most divergent
#      heterozygosity profiles between populations.
#
# Input:  output/LR/LR_<pop>_<rel>_n1000_chunk*.csv
#         (single-locus LR files produced by lr_wrapper.R)
#
# Usage:
#   Rscript code/analyze_locus_inflation.R [output_dir]
#
#   output_dir  Where to write results (default: output/locus_inflation_<date>)
#
# Parallelization: processes each population x relationship combination in
#   parallel using mclapply. Number of cores read from SLURM_CPUS_PER_TASK
#   environment variable (set by --cpus-per-task in SLURM); defaults to 1
#   if run interactively outside SLURM.
#
#   Within each combination, all chunks are read together so true medians
#   are computed across all pairs — not a median of chunk-level medians.
#
# Date: 2026-05-07
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(parallel)
})

log_message <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time()), msg))

# =============================================================================
# ARGUMENT PARSING
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[1] else {
  paste0("output/locus_inflation_", format(Sys.time(), "%Y%m%d_%H%M%S"))
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read cores from SLURM environment; fall back to 1 for interactive use
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))

log_message(paste("Output directory:", output_dir))
log_message(paste("Cores available: ", n_cores))

# =============================================================================
# CONSTANTS
# =============================================================================

POPULATIONS   <- c("AfAm", "Cauc", "Hispanic", "Asian")  # exclude "all"
RELATIONSHIPS <- c("parent_child", "full_siblings", "half_siblings",
                   "cousins", "second_cousins", "unrelated")
LR_DIR        <- "output/LR"

# =============================================================================
# FIND AND GROUP FILES BY POPULATION x RELATIONSHIP
# =============================================================================

log_message("Scanning for single-locus LR files...")

all_lr_files <- list.files(
  LR_DIR,
  pattern = "^LR_(AfAm|Cauc|Hispanic|Asian)_.*_n1000_chunk.*\\.csv$",
  full.names = TRUE
)

if (length(all_lr_files) == 0) {
  stop(paste("No single-locus LR files found in:", LR_DIR,
             "\nExpected pattern: LR_<pop>_<rel>_n1000_chunk*.csv"))
}

log_message(sprintf("Found %d single-locus LR files.", length(all_lr_files)))

# Build a lookup table: file -> pop, rel
file_meta <- data.table(filepath = all_lr_files)
file_meta[, basename          := basename(filepath)]
file_meta[, population        := sub("^LR_(AfAm|Cauc|Hispanic|Asian)_.*", "\\1", basename)]
file_meta[, known_relationship := sub("^LR_(?:AfAm|Cauc|Hispanic|Asian)_(.+)_n1000_chunk.*", "\\1", basename)]

# Group into a named list: one entry per pop x rel combination
combo_groups <- split(file_meta$filepath,
                      paste(file_meta$population,
                            file_meta$known_relationship, sep = "_"))

log_message(sprintf("Processing %d population x relationship combinations across %d cores.",
                    length(combo_groups), n_cores))

# =============================================================================
# PROCESSING FUNCTION — one population x relationship combination
#
# Reads ALL chunks for this combination at once so true medians are computed
# across all pairs, not approximated from chunk-level summaries.
# =============================================================================

process_combination <- function(files, combo_label) {

  cat(sprintf("  [%s] Reading %d chunks...\n", combo_label, length(files)))

  # Read all chunks for this combination and bind
  dt <- rbindlist(lapply(files, fread))

  # Keep only named populations (drop "all" tested population)
  dt <- dt[population %in% c("AfAm", "Cauc", "Hispanic", "Asian") &
           tested_population %in% c("AfAm", "Cauc", "Hispanic", "Asian")]

  # -------------------------------------------------------------------------
  # INFLATION: log10(LR_wrong / LR_correct) per pair x locus
  #
  # Restrict to true positives: tested_relationship == known_relationship.
  # For each pair x locus, join the matched LR (tested_pop == true_pop)
  # against each mismatched LR (tested_pop != true_pop).
  # -------------------------------------------------------------------------

  true_pos <- dt[known_relationship == tested_relationship]

  if (nrow(true_pos) == 0) {
    cat(sprintf("  [%s] No true positives found, skipping inflation.\n", combo_label))
    inflation_summary <- NULL
  } else {

    correct_lr <- true_pos[
      tested_population == population,
      .(pair_id, population, known_relationship, locus, LR_correct = LR)
    ]

    wrong_lr <- true_pos[
      tested_population != population,
      .(pair_id, population, known_relationship, locus, tested_population, LR_wrong = LR)
    ]

    ratio_dt <- merge(wrong_lr, correct_lr,
                      by = c("pair_id", "population", "known_relationship", "locus"))

    # Replace zero/negative LRs with NA before log transform
    ratio_dt[LR_correct <= 0, LR_correct := NA]
    ratio_dt[LR_wrong   <= 0, LR_wrong   := NA]
    ratio_dt[, log10_ratio := log10(LR_wrong) - log10(LR_correct)]

    # True median across all pairs (all chunks read together)
    inflation_summary <- ratio_dt[
      !is.na(log10_ratio),
      .(
        median_log10_ratio = median(log10_ratio),
        mean_log10_ratio   = mean(log10_ratio),
        sd_log10_ratio     = sd(log10_ratio),
        ci_lower           = quantile(log10_ratio, 0.025),
        ci_upper           = quantile(log10_ratio, 0.975),
        n_pairs            = .N
      ),
      by = .(locus, population, tested_population, known_relationship)
    ]

    rm(correct_lr, wrong_lr, ratio_dt, true_pos)
  }

  # -------------------------------------------------------------------------
  # HETEROZYGOSITY: proportion of individuals where focal_A1 != focal_A2
  #
  # Deduplicate to one row per pair x locus first — genotype columns are
  # identical across all tested_relationship x tested_population rows for
  # the same pair x locus.
  # -------------------------------------------------------------------------

  geno <- unique(dt[, .(pair_id, population, locus, focal_A1, focal_A2)])
  geno[, is_het := as.integer(as.character(focal_A1) != as.character(focal_A2))]

  heteroz_summary <- geno[,
    .(
      obs_heterozygosity = mean(is_het, na.rm = TRUE),
      n_individuals      = .N
    ),
    by = .(locus, population)
  ]

  rm(geno, dt)

  list(inflation = inflation_summary, heterozygosity = heteroz_summary)
}

# =============================================================================
# RUN IN PARALLEL ACROSS COMBINATIONS
# =============================================================================

log_message("Starting parallel processing...")

results <- mclapply(
  X              = names(combo_groups),
  FUN            = function(combo_label) {
    process_combination(files       = combo_groups[[combo_label]],
                        combo_label = combo_label)
  },
  mc.cores       = n_cores,
  mc.preschedule = FALSE   # dynamic scheduling since combinations vary in file count
)

names(results) <- names(combo_groups)

log_message("Parallel processing complete. Combining results...")

# =============================================================================
# COMBINE AND WRITE OUTPUT 1: LOCUS INFLATION SUMMARY
# =============================================================================

inflation_list <- lapply(results, `[[`, "inflation")
inflation_all  <- rbindlist(Filter(Negate(is.null), inflation_list))

# Rank loci within each true_pop x tested_pop x relationship group
# (rank 1 = locus contributing most inflation)
inflation_all[,
  locus_rank := frank(-median_log10_ratio, ties.method = "min"),
  by = .(population, tested_population, known_relationship)
]

setorder(inflation_all, population, tested_population, known_relationship, locus_rank)

inflation_file <- file.path(output_dir, "locus_inflation_summary.csv")
fwrite(inflation_all, inflation_file)
log_message(sprintf("Inflation summary written: %s (%d rows)", inflation_file, nrow(inflation_all)))

# =============================================================================
# COMBINE AND WRITE OUTPUT 2: LOCUS HETEROZYGOSITY SUMMARY
# =============================================================================

heteroz_list <- lapply(results, `[[`, "heterozygosity")
heteroz_all  <- rbindlist(Filter(Negate(is.null), heteroz_list))

# Pool across relationship combinations — same locus x population appears
# in each relationship's results; weight by n_individuals when pooling
heteroz_pooled <- heteroz_all[,
  .(
    obs_heterozygosity = weighted.mean(obs_heterozygosity, n_individuals),
    n_individuals      = sum(n_individuals)
  ),
  by = .(locus, population)
]

# Wide format: one row per locus, one column per population — easier to scan
heteroz_wide <- dcast(heteroz_pooled, locus ~ population,
                      value.var = "obs_heterozygosity")

# Add AfAm - Asian difference (the key comparison)
if (all(c("AfAm", "Asian") %in% names(heteroz_wide))) {
  heteroz_wide[, AfAm_minus_Asian := AfAm - Asian]
  setorder(heteroz_wide, -abs(AfAm_minus_Asian))
}

heteroz_file <- file.path(output_dir, "locus_heterozygosity_summary.csv")
fwrite(heteroz_wide, heteroz_file)
log_message(sprintf("Heterozygosity summary written: %s (%d rows)", heteroz_file, nrow(heteroz_wide)))

# =============================================================================
# CONSOLE SUMMARY
# =============================================================================

cat("\n")
cat("=============================================================================\n")
cat("TOP 5 LOCI BY INFLATION — AfAm true pairs tested with Asian frequencies\n")
cat("  (full_siblings)\n")
cat("=============================================================================\n")

top5 <- inflation_all[
  population == "AfAm" & tested_population == "Asian" &
  known_relationship == "full_siblings"
][order(locus_rank)][1:min(5, .N)]

print(top5[, .(locus, locus_rank, median_log10_ratio, ci_lower, ci_upper, n_pairs)])

cat("\n")
cat("=============================================================================\n")
cat("HETEROZYGOSITY — AfAm vs Asian (largest absolute difference first)\n")
cat("=============================================================================\n")

if (all(c("AfAm", "Asian") %in% names(heteroz_wide))) {
  print(head(heteroz_wide[, .(locus, AfAm, Asian, AfAm_minus_Asian)], 10))
}

cat("\n")
cat("OUTPUT FILES:\n")
cat(sprintf("  %s\n", inflation_file))
cat(sprintf("  %s\n", heteroz_file))
cat("=============================================================================\n")

log_message("Done.")
