# ------------------------------------------------------------------------------
# Module 4 (fast path): Vectorized Single Locus LR Calculator for the ranking test
# ------------------------------------------------------------------------------
#
# Why this exists:
#   calculate_single_locus_lr() (module4_single_locus_LR.R) loops row-by-row
#   over pair_data and calls kinship_calculation() once per row. Inside that
#   call, both kinship_coefficients and the 14k-row allele_frequency_data get
#   re-filtered with dplyr::filter() EVERY row, plus a tibble is rebuilt via
#   tidy-eval (!!!row) every row. For the main pipeline (~1000 pairs/file)
#   that overhead is tolerable. For the focal ranking test, one replicate
#   fans out to 10,000 unrelated-pool candidates x ~13-20 loci x 2 tested
#   relationships = 260k-400k of these calls PER PAIR, which is what turned
#   a single pair into a 44-minute job.
#
#   This file does not touch module4_single_locus_LR.R or
#   LR_kinship_utility_functions.R, so the main pipeline (lr_wrapper.R,
#   combined_lr_wrapper.R, and the already-generated output/LR and
#   output/combined_LR files) is completely unaffected. This is only wired
#   into the ranking test path.
#
# Validated against the original row-by-row logic (shared-allele counting +
# LR formula) across all 2401 possible 2-allele genotype combinations from
# a 7-allele pool: 0 mismatches. See conversation for the check script.
# Still recommend spot-checking against a known-good run before trusting
# at scale - see validation note at the bottom of this file.
#
# ------------------------------------------------------------------------------

library(data.table)

FALLBACK_FREQ_FAST <- 5 / (2 * 1036)   # must match FALLBACK_FREQ in LR_kinship_utility_functions.R

#' Vectorized single-locus LR calculation (drop-in for calculate_single_locus_lr)
#'
#' Same required input columns and same required output columns as
#' calculate_single_locus_lr() (module4_single_locus_LR.R), so it can be
#' substituted directly inside calculate_ranking_lrs() (module 11).
#' Omits genotype_match (not used by calculate_combined_lr / module 5).
#'
#' @param pair_data Data frame/data.table with columns: batch_id, pair_id,
#'        population, locus, focal_A1, focal_A2, ind2_A1, ind2_A2,
#'        known_relationship
#' @param tested_relationship Character scalar, e.g. "full_siblings"
#' @param tested_populations Character vector, populations to test LRs against
#' @param allele_frequency_data Data frame with population, marker, allele, frequency
#' @param kinship_coefficients Data frame with relationship_type, k0, k1, k2
#' @return data.table with one row per pair_id x locus x tested_population,
#'         columns: batch_id, pair_id, population, known_relationship, locus,
#'         focal_A1, focal_A2, ind2_A1, ind2_A2, shared_alleles,
#'         tested_relationship, tested_population, LR
calculate_single_locus_lr_fast <- function(pair_data,
                                            tested_relationship,
                                            tested_populations,
                                            allele_frequency_data,
                                            kinship_coefficients) {

  required_cols <- c("batch_id", "pair_id", "population", "locus",
                      "focal_A1", "focal_A2", "ind2_A1", "ind2_A2",
                      "known_relationship")
  missing_cols <- setdiff(required_cols, names(pair_data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }

  if (!tested_relationship %in% kinship_coefficients$relationship_type) {
    stop(paste("tested_relationship", tested_relationship,
               "not found in kinship_coefficients"))
  }

  # Kinship coefficients looked up ONCE per call (was re-filtered every row)
  kin_row <- kinship_coefficients[kinship_coefficients$relationship_type == tested_relationship, ]
  if (nrow(kin_row) != 1L) {
    stop("tested_relationship not found or duplicated in kinship_coefficients")
  }
  k0 <- kin_row$k0[1]; k1 <- kin_row$k1[1]; k2 <- kin_row$k2[1]

  dt <- as.data.table(pair_data)
  dt[, `:=`(
    A1 = as.character(focal_A1),
    A2 = as.character(focal_A2),
    B1 = as.character(ind2_A1),
    B2 = as.character(ind2_A2)
  )]

  # ---------------------------------------------------------------------
  # Shared alleles (with multiplicity), fully vectorized
  # ---------------------------------------------------------------------
  dt[, ind1_hom := A1 == A2]
  dt[, ind2_hom := B1 == B2]

  dt[, shared_alleles := fcase(
    ind1_hom & ind2_hom,   as.integer(A1 == B1) * 2L,
    ind1_hom & !ind2_hom,  as.integer(A1 == B1 | A1 == B2),
    !ind1_hom & ind2_hom,  as.integer(B1 == A1 | B1 == A2),
    default =              as.integer(A1 == B1 | A1 == B2) + as.integer(A2 == B1 | A2 == B2)
  )]

  # The specific shared allele value, only meaningful when shared_alleles == 1
  # NOTE: must be a row-wise comparison (A1 == B1 | A1 == B2), NOT
  # `A1 %chin% c(B1, B2)` - c(B1, B2) concatenates the whole columns
  # together rather than pairing B1[i]/B2[i] for row i, which silently
  # matches against every OTHER row's alleles too.
  dt[, shared_val := fifelse(A1 == B1 | A1 == B2, A1, A2)]

  # Rxp multiplier for the 1-shared-allele case: 2 if either individual is
  # homozygous ("AA-AB"/"AB-AA"), 4 if both heterozygous ("AB-AC")
  dt[, mult := fifelse(ind1_hom | ind2_hom, 2, 4)]

  # ---------------------------------------------------------------------
  # Expand to one row per tested_population (was per-row map_dfr before)
  # ---------------------------------------------------------------------
  n_pop <- length(tested_populations)
  dt <- dt[rep(seq_len(.N), each = n_pop)]
  dt[, tested_population := rep(tested_populations, times = nrow(pair_data))]

  # ---------------------------------------------------------------------
  # Single hash join for allele frequencies (was a dplyr::filter() of the
  # full 14k-row allele_frequency_data table on EVERY row before)
  # ---------------------------------------------------------------------
  freq_key <- as.data.table(allele_frequency_data)[, .(population, marker, allele, frequency)]
  setnames(freq_key, c("population", "marker", "allele"),
                     c("tested_population", "locus", "allele_lookup"))

  # which (tested_population, locus) combos exist at all, to distinguish
  # "locus missing entirely -> LR NA" from "allele missing -> contingency freq"
  combo_exists <- unique(freq_key[, .(tested_population, locus)])
  combo_exists[, locus_present := TRUE]

  dt[, allele_A_lookup := fifelse(shared_alleles == 1, shared_val, A1)]
  dt[, allele_B_lookup := A2]  # only used when shared_alleles == 2 & heterozygous

  dt <- merge(dt, combo_exists, by = c("tested_population", "locus"), all.x = TRUE, sort = FALSE)
  dt[is.na(locus_present), locus_present := FALSE]

  dt <- merge(dt, freq_key,
              by.x = c("tested_population", "locus", "allele_A_lookup"),
              by.y = c("tested_population", "locus", "allele_lookup"),
              all.x = TRUE, sort = FALSE)
  setnames(dt, "frequency", "pA")

  dt <- merge(dt, freq_key,
              by.x = c("tested_population", "locus", "allele_B_lookup"),
              by.y = c("tested_population", "locus", "allele_lookup"),
              all.x = TRUE, sort = FALSE)
  setnames(dt, "frequency", "pB")

  # Contingency fallback (mirrors fetch_freq() in LR_kinship_utility_functions.R):
  # only apply the fallback frequency when the locus/population combo exists
  # but this particular allele wasn't observed. If the whole locus is missing
  # for this population, leave NA so LR comes out NA downstream.
  dt[locus_present == TRUE & (is.na(pA) | pA == 0), pA := FALLBACK_FREQ_FAST]
  dt[locus_present == TRUE & shared_alleles == 2 & !ind1_hom & (is.na(pB) | pB == 0),
     pB := FALLBACK_FREQ_FAST]

  # ---------------------------------------------------------------------
  # Vectorized LR formula (mirrors calculate_likelihood_ratio())
  # ---------------------------------------------------------------------
  dt[, LR := fcase(
    shared_alleles == 0,                 k0,
    shared_alleles == 1,                 k0 + (k1 / (mult * pA)),
    shared_alleles == 2 & ind1_hom,      k0 + (k1 / pA) + (k2 / (pA^2)),
    shared_alleles == 2 & !ind1_hom,     k0 + (k1 / ((4 * pA * pB) / (pA + pB))) + (k2 / (2 * pA * pB))
  )]

  dt[, tested_relationship := tested_relationship]

  final_results <- dt[, .(batch_id, pair_id, population, known_relationship, locus,
                           focal_A1, focal_A2, ind2_A1, ind2_A2,
                           shared_alleles, tested_relationship, tested_population, LR)]

  return(final_results)
}

# ------------------------------------------------------------------------------
# VALIDATION NOTE:
#
#   v1 of this file had a real bug: the "both individuals heterozygous"
#   branch of shared_alleles, and the shared_val line, used
#   `A1 %in% c(B1, B2)` / `A1 %chin% c(B1, B2)`. When B1/B2 are whole
#   data.table columns (not scalars), c(B1, B2) concatenates the ENTIRE
#   columns rather than pairing B1[i]/B2[i] for row i - so it was checking
#   whether A1 matched ANY candidate's allele anywhere in the pool, not just
#   the current row's candidate. Caught by diffing against
#   calculate_single_locus_lr() on a 20-candidate slice of real data
#   (compare_old_vs_fast_lr.R): 223/400 rows mismatched, all in the
#   both-heterozygous branch. Fixed by using row-wise `A1 == B1 | A1 == B2`
#   comparisons instead.
#
#   The underlying shared-allele-counting/LR-formula MATH (elementwise
#   comparisons, as now used throughout) was cross-checked in Python against
#   a direct translation of count_shared_alleles() + label_and_genotype() +
#   calculate_likelihood_ratio() across all 2401 possible 2-allele genotype
#   combinations from a 7-allele pool: 0 mismatches. That was never the
#   problem - the problem was an R vectorization footgun introduced when
#   translating the validated math into data.table column operations.
#
#   Before trusting this at scale again: re-run compare_old_vs_fast_lr.R
#   (or equivalent) on real data and confirm 0 mismatches before pointing
#   a real array job at it. Do not assume a second bug-free result on faith -
#   the first "0 mismatches" claim in this file was about formula logic in
#   isolation, not this R implementation, and that distinction is exactly
#   what got missed.
#
# ------------------------------------------------------------------------------
