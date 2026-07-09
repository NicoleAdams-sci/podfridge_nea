# ------------------------------------------------------------------------------
# Module 12: Ranking and Outcome Recorder
# ------------------------------------------------------------------------------
#
# Purpose:
#   Takes combined LR results from Module 11, ranks all database members by
#   combined LR for each focal individual replicate, and records whether the
#   true relative falls within the top N (default top 200).
#
#   Summarizes outcomes across all replicates to answer:
#     - How often does the true relative appear in the top 200?
#     - What is the distribution of true relative ranks?
#     - How does this vary by relationship type, loci set, and population frequency?
#
# Dependencies:
#   - Module 11 output: combined LRs across all replicates
#
# Output:
#   - output/focal_ranking_test/ranking_outcomes_<datetime>.csv
#       One row per replicate x loci_set x tested_relationship x tested_population
#   - output/focal_ranking_test/ranking_summary_<datetime>.csv
#       Aggregated summary across replicates
#
# ------------------------------------------------------------------------------

library(dplyr)
library(data.table)


#' Rank database members and record true relative outcome for all replicates
#'
#' For each replicate (focal individual), loci set, tested relationship, and
#' tested population: ranks all database members by combined LR (descending),
#' finds the true relative's rank, and records whether it falls within top_n.
#'
#' @param lr_results Data frame, output from Module 11
#'        calculate_ranking_lrs_replicates(). Must contain columns:
#'        focal_id, pair_id, individual_id, loci_set, tested_relationship,
#'        tested_population, known_relationship, combined_LR, is_true_relative.
#' @param top_n Integer, rank threshold to evaluate (default 200).
#' @param output_dir Character, directory to save results.
#' @param save_results Logical, whether to save results to CSV (default TRUE).
#' @return List with two data frames:
#'         $outcomes  — one row per replicate x loci_set x tested_relationship
#'                      x tested_population, with rank, tied_group_size, and
#'                      in_top_n
#'         $summary   — aggregated proportion in top_n across replicates

record_ranking_outcomes <- function(lr_results,
                                    top_n = 200,
                                    output_dir = "output/focal_ranking_test",
                                    save_results = TRUE) {

  # --- Input validation ---
  required_cols <- c("focal_id", "pair_id", "individual_id", "loci_set",
                     "tested_relationship", "tested_population",
                     "known_relationship", "combined_LR", "is_true_relative")
  missing_cols <- setdiff(required_cols, names(lr_results))
  if (length(missing_cols) > 0)
    stop(paste("lr_results missing columns:", paste(missing_cols, collapse = ", ")))

  n_replicates <- length(unique(lr_results$focal_id))
  cat(sprintf("Recording ranking outcomes for %d replicates (top_%d threshold)...\n",
              n_replicates, top_n))

  # --- Rank within each replicate x loci_set x tested_relationship x tested_population ---
  # For each grouping: rank all pairs by combined_LR descending,
  # then find the true relative's rank.
  #
  # Note: ties broken by average rank (standard competition ranking).
  # The true relative row has is_true_relative == TRUE.

  outcomes <- lr_results |>
    group_by(focal_id, loci_set, tested_relationship, tested_population) |>
    mutate(
      rank = rank(-combined_LR, ties.method = "average"),
      n_database = n_distinct(individual_id),
      # How many candidates (unrelateds + the true relative itself) landed
      # at EXACTLY the same combined_LR as the true relative in this group.
      # A large tied_group_size (common when combined_LR == 0, e.g. testing
      # a full/half-sib true relative against the parent_child hypothesis,
      # which has k0 = 0) means "rank" is an average over a big block of
      # indistinguishable candidates, not a real, singular position.
      tied_group_size = sum(combined_LR == combined_LR[is_true_relative][1], na.rm = TRUE)
    ) |>
    ungroup() |>
    filter(is_true_relative == TRUE) |>
    mutate(
      in_top_n   = rank <= top_n,
      top_n_used = top_n
    ) |>
    select(focal_id, loci_set, tested_relationship, tested_population,
           known_relationship, individual_id, combined_LR,
           rank, n_database, tied_group_size, in_top_n, top_n_used)

  # --- Summarize across replicates ---
  summary_results <- outcomes |>
    group_by(loci_set, tested_relationship, tested_population, known_relationship) |>
    summarize(
      n_replicates      = n(),
      n_in_top_n        = sum(in_top_n),
      prop_in_top_n     = mean(in_top_n),
      mean_rank         = mean(rank),
      median_rank       = median(rank),
      sd_rank           = sd(rank),
      min_rank          = min(rank),
      max_rank          = max(rank),
      mean_tied_group_size   = mean(tied_group_size),
      median_tied_group_size = median(tied_group_size),
      top_n_used        = first(top_n_used),
      .groups = "drop"
    )

  # Print summary to console
  cat("\n--- Ranking Summary ---\n")
  print(as.data.frame(summary_results |>
    select(loci_set, tested_relationship, tested_population,
           n_replicates, prop_in_top_n, median_rank)))

  # --- Save results ---
  if (save_results) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    datetime_str <- format(Sys.time(), "%Y%m%d_%H%M%S")

    outcomes_path <- file.path(output_dir,
                               paste0("ranking_outcomes_", datetime_str, ".csv"))
    summary_path  <- file.path(output_dir,
                               paste0("ranking_summary_", datetime_str, ".csv"))

    fwrite(outcomes, outcomes_path)
    fwrite(summary_results, summary_path)

    cat(sprintf("\nOutcomes saved:  %s\n", outcomes_path))
    cat(sprintf("Summary saved:   %s\n", summary_path))
  }

  return(list(
    outcomes = as.data.frame(outcomes),
    summary  = as.data.frame(summary_results)
  ))
}


# ------------------------------------------------------------------------------
# Usage Examples (commented out)
# ------------------------------------------------------------------------------
#
# # --- Run ranking and record outcomes ---
# results <- record_ranking_outcomes(
#   lr_results  = all_lr_results,   # from Module 11
#   top_n       = 200,
#   output_dir  = "output/focal_ranking_test"
# )
#
# # --- Access outcomes and summary ---
# outcomes_df <- results$outcomes   # one row per replicate x condition
# summary_df  <- results$summary    # aggregated proportions across replicates
#
# # --- Check a quick result ---
# # What proportion of parent_child relatives ranked in top 200 using core_13?
# summary_df |>
#   filter(loci_set == "core_13", tested_relationship == "parent_child")
