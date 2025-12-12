# Test Module 9: Statistics functions for combined LR analysis
# - calculate_summary_stats()
# - calculate_ratio_stats()
# - calculate_cutoffs()
# - calculate_proportions_exceeding_cutoffs()

library(testthat)
library(dplyr)
library(data.table)

# =============================================================================
# Define functions locally (copied from module9_combinedLR_stats_functions.R)
# =============================================================================

calculate_summary_stats <- function(combined_lrs) {
  summary_stats <- combined_lrs |>
    group_by(known_relationship, population, loci_set, tested_population, tested_relationship, is_correct_pop) |>
    summarize(
      n = n(),
      mean_LR = mean(combined_LR, na.rm = TRUE),
      median_LR = median(combined_LR, na.rm = TRUE),
      sd_LR = sd(combined_LR, na.rm = TRUE),
      min_LR = min(combined_LR, na.rm = TRUE),
      max_LR = max(combined_LR, na.rm = TRUE),
      lower_95 = quantile(combined_LR, 0.025, na.rm = TRUE),
      upper_95 = quantile(combined_LR, 0.975, na.rm = TRUE),
      .groups = 'drop'
    ) |>
    ungroup()
  return(summary_stats)
}

calculate_cutoffs <- function(all_combined, fp_rates) {
  cutoffs <- all_combined %>%
    filter(known_relationship == "unrelated", is_correct_pop == TRUE) %>% 
    group_by(loci_set, tested_population, tested_relationship) %>%
    summarize(
      fixed_cutoff = 1.00,
      cutoff_1 = quantile(combined_LR, probs = 1 - fp_rates[1] / 100, na.rm = TRUE),
      cutoff_0_1 = quantile(combined_LR, probs = 1 - fp_rates[2] / 100, na.rm = TRUE),
      cutoff_0_01 = quantile(combined_LR, probs = 1 - fp_rates[3] / 100, na.rm = TRUE),
      n_unrelated = n(),
      .groups = 'drop'
    )
  return(cutoffs)
}

calculate_proportions_exceeding_cutoffs <- function(all_combined, cutoffs) {
  df_with_cutoffs <- left_join(
    all_combined,
    cutoffs,
    by = c("loci_set", "tested_population", "tested_relationship")
  )
  df_with_cutoffs <- df_with_cutoffs %>%
    mutate(
      exceeds_fixed_cutoff = combined_LR > fixed_cutoff,
      exceeds_cutoff_1     = combined_LR > cutoff_1,
      exceeds_cutoff_0_1   = combined_LR > cutoff_0_1,
      exceeds_cutoff_0_01  = combined_LR > cutoff_0_01
    )
  proportions_exceeding <- df_with_cutoffs %>%
    group_by(
      population, known_relationship, tested_relationship,
      loci_set, tested_population, is_correct_pop
    ) %>%
    summarize(
      proportion_exceeding_fixed = sum(exceeds_fixed_cutoff, na.rm = TRUE) / n(),
      proportion_exceeding_1     = sum(exceeds_cutoff_1, na.rm = TRUE) / n(),
      proportion_exceeding_0_1   = sum(exceeds_cutoff_0_1, na.rm = TRUE) / n(),
      proportion_exceeding_0_01  = sum(exceeds_cutoff_0_01, na.rm = TRUE) / n(),
      n_related = n(),
      .groups = 'drop'
    )
  return(proportions_exceeding)
}

calculate_ratio_stats <- function(all_combined) {
  all_combined <- as.data.table(all_combined)
  combined_lrs_correct <- all_combined[
    is_correct_pop == TRUE & (tested_population == population),
    .(batch_id, pair_id, population, known_relationship, tested_relationship, loci_set, correct_LR = combined_LR)
  ]
  combined_lrs_wrong <- all_combined[
    is_correct_pop == FALSE & (tested_population != population),
    .(batch_id, pair_id, population, known_relationship, tested_relationship, loci_set, tested_population, wrong_LR = combined_LR)
  ]
  combined_lrs_correct_unique <- unique(combined_lrs_correct,
                                        by = c("batch_id", "pair_id", "population",
                                               "known_relationship", "tested_relationship", "loci_set"))
  combined_lrs_ratio <- merge(
    combined_lrs_wrong,
    combined_lrs_correct_unique,
    by = c("batch_id", "pair_id", "population", "known_relationship", "loci_set", "tested_relationship")
  )
  combined_lrs_ratio[, ratio := wrong_LR / correct_LR]
  ratio_summary <- combined_lrs_ratio %>%
    group_by(population, known_relationship, tested_relationship, loci_set, tested_population) %>%
    summarize(
      n = n(),
      mean_ratio = mean(ratio, na.rm = TRUE),
      median_ratio = median(ratio, na.rm = TRUE),
      sd_ratio = sd(ratio, na.rm = TRUE),
      min_ratio = min(ratio, na.rm = TRUE),
      max_ratio = max(ratio, na.rm = TRUE),
      lower_95 = quantile(ratio, 0.025, na.rm = TRUE),
      upper_95 = quantile(ratio, 0.975, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    ungroup()
  return(list(ratio_summary = as.data.frame(ratio_summary),
              combined_lrs_ratio = as.data.frame(combined_lrs_ratio)))
}

# =============================================================================
# MOCK DATA SETUP
# =============================================================================

# Create mock combined LR data with known values for testing
# 10 pairs: 5 parent_child, 5 unrelated
# Each tested with parent_child hypothesis
create_mock_combined_data <- function() {
  data.frame(
    batch_id = rep("batch1", 20),
    pair_id = c(rep(paste0("pc_", 1:5), 2), rep(paste0("un_", 1:5), 2)),
    population = rep("PopA", 20),
    known_relationship = c(rep("parent_child", 10), rep("unrelated", 10)),
    loci_set = rep("set1", 20),
    tested_relationship = rep("parent_child", 20),
    tested_population = c(rep("PopA", 5), rep("PopB", 5), rep("PopA", 5), rep("PopB", 5)),
    is_correct_pop = c(rep(TRUE, 5), rep(FALSE, 5), rep(TRUE, 5), rep(FALSE, 5)),
    # Parent-child pairs: high LRs when tested correctly
    # Unrelated pairs: LRs around 1
    combined_LR = c(
      # parent_child pairs, tested with PopA (correct): high LRs
      100, 150, 200, 120, 130,
      # parent_child pairs, tested with PopB (wrong): slightly different
      90, 140, 180, 110, 125,
      # unrelated pairs, tested with PopA (correct): around 1
      0.8, 1.0, 1.2, 0.9, 1.1,
      # unrelated pairs, tested with PopB (wrong): slightly different
      0.7, 0.95, 1.3, 0.85, 1.15
    ),
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# TEST: calculate_summary_stats()
# =============================================================================

test_that("calculate_summary_stats returns correct mean for known data", {
  mock_data <- create_mock_combined_data()
  result <- calculate_summary_stats(mock_data)
  
  # Parent-child pairs with correct population
  pc_correct <- result |> 
    filter(known_relationship == "parent_child", is_correct_pop == TRUE)
  
  # Mean of 100, 150, 200, 120, 130 = 700/5 = 140
  expect_equal(pc_correct$mean_LR, 140)
  expect_equal(pc_correct$n, 5)
})

test_that("calculate_summary_stats returns correct median", {
  mock_data <- create_mock_combined_data()
  result <- calculate_summary_stats(mock_data)
  
  # Parent-child, correct pop: sorted = 100, 120, 130, 150, 200 -> median = 130
  pc_correct <- result |> 
    filter(known_relationship == "parent_child", is_correct_pop == TRUE)
  expect_equal(pc_correct$median_LR, 130)
})

test_that("calculate_summary_stats returns correct min and max", {
  mock_data <- create_mock_combined_data()
  result <- calculate_summary_stats(mock_data)
  
  pc_correct <- result |> 
    filter(known_relationship == "parent_child", is_correct_pop == TRUE)
  
  expect_equal(pc_correct$min_LR, 100)
  expect_equal(pc_correct$max_LR, 200)
})

test_that("calculate_summary_stats groups correctly", {
  mock_data <- create_mock_combined_data()
  result <- calculate_summary_stats(mock_data)
  
  # Should have 4 rows: 2 relationships × 2 is_correct_pop values
  expect_equal(nrow(result), 4)
})

test_that("calculate_summary_stats handles NA values", {
  mock_data <- create_mock_combined_data()
  mock_data$combined_LR[1] <- NA
  
  result <- calculate_summary_stats(mock_data)
  
  # Should still compute stats (na.rm = TRUE)
  pc_correct <- result |> 
    filter(known_relationship == "parent_child", is_correct_pop == TRUE)
  
  # Mean of 150, 200, 120, 130 = 600/4 = 150 (excluding NA)
  expect_equal(pc_correct$mean_LR, 150)
  expect_equal(pc_correct$n, 5)  # n() still counts all rows
})

# =============================================================================
# TEST: calculate_cutoffs()
# =============================================================================

test_that("calculate_cutoffs uses only unrelated pairs with correct pop", {
  mock_data <- create_mock_combined_data()
  fp_rates <- c(1, 0.1, 0.01)
  
  result <- calculate_cutoffs(mock_data, fp_rates)
  
  # Should only use the 5 unrelated pairs with is_correct_pop == TRUE
  expect_equal(result$n_unrelated, 5)
})

test_that("calculate_cutoffs returns fixed_cutoff of 1", {
  mock_data <- create_mock_combined_data()
  fp_rates <- c(1, 0.1, 0.01)
  
  result <- calculate_cutoffs(mock_data, fp_rates)
  
  expect_equal(result$fixed_cutoff, 1.0)
})

test_that("calculate_cutoffs calculates correct quantiles", {
  # Create data with known quantiles
  mock_simple <- data.frame(
    batch_id = rep("b1", 100),
    pair_id = paste0("p", 1:100),
    population = rep("PopA", 100),
    known_relationship = rep("unrelated", 100),
    loci_set = rep("set1", 100),
    tested_relationship = rep("parent_child", 100),
    tested_population = rep("PopA", 100),
    is_correct_pop = rep(TRUE, 100),
    # LRs from 1 to 100
    combined_LR = 1:100,
    stringsAsFactors = FALSE
  )
  
  fp_rates <- c(1, 0.1, 0.01)
  result <- calculate_cutoffs(mock_simple, fp_rates)
  
  # 1% FPR = 99th percentile
  # For 1:100, quantile(0.99) ≈ 99.01
  expect_equal(result$cutoff_1, quantile(1:100, 0.99), tolerance = 0.1)
})

test_that("calculate_cutoffs groups by loci_set and tested_population", {
  mock_data <- create_mock_combined_data()
  
  # Add another loci_set
  mock_data2 <- mock_data
  mock_data2$loci_set <- "set2"
  mock_multi <- rbind(mock_data, mock_data2)
  
  fp_rates <- c(1, 0.1, 0.01)
  result <- calculate_cutoffs(mock_multi, fp_rates)
  
  # Should have 2 rows (one per loci_set, same tested_population)
  expect_equal(nrow(result), 2)
})

# =============================================================================
# TEST: calculate_proportions_exceeding_cutoffs()
# =============================================================================

test_that("calculate_proportions_exceeding_cutoffs correctly identifies exceedances", {
  # Simple case: all LRs exceed fixed cutoff of 1.0
  mock_data <- data.frame(
    batch_id = rep("b1", 5),
    pair_id = paste0("p", 1:5),
    population = rep("PopA", 5),
    known_relationship = rep("parent_child", 5),
    loci_set = rep("set1", 5),
    tested_relationship = rep("parent_child", 5),
    tested_population = rep("PopA", 5),
    is_correct_pop = rep(TRUE, 5),
    combined_LR = c(10, 20, 30, 40, 50),  # All > 1
    stringsAsFactors = FALSE
  )
  
  cutoffs <- data.frame(
    loci_set = "set1",
    tested_population = "PopA",
    tested_relationship = "parent_child",
    fixed_cutoff = 1.0,
    cutoff_1 = 25,      # 3/5 exceed this (30, 40, 50)
    cutoff_0_1 = 45,    # 1/5 exceeds this (50)
    cutoff_0_01 = 55,   # 0/5 exceed this
    stringsAsFactors = FALSE
  )
  
  result <- calculate_proportions_exceeding_cutoffs(mock_data, cutoffs)
  
  expect_equal(result$proportion_exceeding_fixed, 1.0)  # 5/5
  expect_equal(result$proportion_exceeding_1, 0.6)      # 3/5
  expect_equal(result$proportion_exceeding_0_1, 0.2)    # 1/5
  expect_equal(result$proportion_exceeding_0_01, 0.0)   # 0/5
  expect_equal(result$n_related, 5)
})

test_that("calculate_proportions_exceeding_cutoffs handles no exceedances", {
  mock_data <- data.frame(
    batch_id = rep("b1", 5),
    pair_id = paste0("p", 1:5),
    population = rep("PopA", 5),
    known_relationship = rep("unrelated", 5),
    loci_set = rep("set1", 5),
    tested_relationship = rep("parent_child", 5),
    tested_population = rep("PopA", 5),
    is_correct_pop = rep(TRUE, 5),
    combined_LR = c(0.5, 0.6, 0.7, 0.8, 0.9),  # All < 1
    stringsAsFactors = FALSE
  )
  
  cutoffs <- data.frame(
    loci_set = "set1",
    tested_population = "PopA",
    tested_relationship = "parent_child",
    fixed_cutoff = 1.0,
    cutoff_1 = 10,
    cutoff_0_1 = 50,
    cutoff_0_01 = 100,
    stringsAsFactors = FALSE
  )
  
  result <- calculate_proportions_exceeding_cutoffs(mock_data, cutoffs)
  
  expect_equal(result$proportion_exceeding_fixed, 0.0)
  expect_equal(result$proportion_exceeding_1, 0.0)
})

test_that("calculate_proportions_exceeding_cutoffs groups correctly", {
  mock_data <- data.frame(
    batch_id = rep("b1", 10),
    pair_id = paste0("p", 1:10),
    population = rep("PopA", 10),
    known_relationship = c(rep("parent_child", 5), rep("unrelated", 5)),
    loci_set = rep("set1", 10),
    tested_relationship = rep("parent_child", 10),
    tested_population = rep("PopA", 10),
    is_correct_pop = rep(TRUE, 10),
    combined_LR = c(100, 200, 300, 400, 500, 0.5, 0.6, 0.7, 0.8, 0.9),
    stringsAsFactors = FALSE
  )
  
  cutoffs <- data.frame(
    loci_set = "set1",
    tested_population = "PopA",
    tested_relationship = "parent_child",
    fixed_cutoff = 1.0,
    cutoff_1 = 10,
    cutoff_0_1 = 50,
    cutoff_0_01 = 100,
    stringsAsFactors = FALSE
  )
  
  result <- calculate_proportions_exceeding_cutoffs(mock_data, cutoffs)
  
  # Should have 2 rows (parent_child and unrelated)
  expect_equal(nrow(result), 2)
  
  # Parent-child: all 5 exceed all cutoffs
  pc_result <- result |> filter(known_relationship == "parent_child")
  expect_equal(pc_result$proportion_exceeding_fixed, 1.0)
  
  # Unrelated: none exceed fixed cutoff
  un_result <- result |> filter(known_relationship == "unrelated")
  expect_equal(un_result$proportion_exceeding_fixed, 0.0)
})

# =============================================================================
# TEST: calculate_ratio_stats()
# =============================================================================

test_that("calculate_ratio_stats computes correct ratio", {
  # Create data with one pair that has both correct and wrong pop results
  mock_data <- data.frame(
    batch_id = rep("b1", 4),
    pair_id = rep("pair1", 4),
    population = rep("PopA", 4),
    known_relationship = rep("parent_child", 4),
    loci_set = rep("set1", 4),
    tested_relationship = rep("parent_child", 4),
    tested_population = c("PopA", "PopB", "PopA", "PopB"),
    is_correct_pop = c(TRUE, FALSE, TRUE, FALSE),
    # Correct (PopA): LR = 100
    # Wrong (PopB): LR = 80
    # Ratio should be 80/100 = 0.8
    combined_LR = c(100, 80, 100, 80),
    stringsAsFactors = FALSE
  )
  
  result <- calculate_ratio_stats(mock_data)
  
  expect_equal(nrow(result$ratio_summary), 1)
  expect_equal(result$ratio_summary$mean_ratio, 0.8)
  expect_equal(result$ratio_summary$median_ratio, 0.8)
})

test_that("calculate_ratio_stats handles multiple pairs", {
  mock_data <- data.frame(
    batch_id = rep("b1", 4),
    pair_id = c("pair1", "pair1", "pair2", "pair2"),
    population = rep("PopA", 4),
    known_relationship = rep("parent_child", 4),
    loci_set = rep("set1", 4),
    tested_relationship = rep("parent_child", 4),
    tested_population = c("PopA", "PopB", "PopA", "PopB"),
    is_correct_pop = c(TRUE, FALSE, TRUE, FALSE),
    combined_LR = c(100, 80,    # pair1: ratio = 0.8
                    200, 100),  # pair2: ratio = 0.5
    stringsAsFactors = FALSE
  )
  
  result <- calculate_ratio_stats(mock_data)
  
  # Mean ratio = (0.8 + 0.5) / 2 = 0.65
  expect_equal(result$ratio_summary$mean_ratio, 0.65)
  expect_equal(result$ratio_summary$n, 2)
  
  # Check raw ratios
  expect_equal(nrow(result$combined_lrs_ratio), 2)
  ratios <- sort(result$combined_lrs_ratio$ratio)
  expect_equal(ratios, c(0.5, 0.8))
})

test_that("calculate_ratio_stats returns both summary and raw data", {
  mock_data <- create_mock_combined_data()
  result <- calculate_ratio_stats(mock_data)
  
  expect_true("ratio_summary" %in% names(result))
  expect_true("combined_lrs_ratio" %in% names(result))
  expect_true(is.data.frame(result$ratio_summary))
  expect_true(is.data.frame(result$combined_lrs_ratio))
})

test_that("calculate_ratio_stats only compares correct vs wrong pop for same pair", {
  # Data with correct pop only - should produce no ratios
  mock_correct_only <- data.frame(
    batch_id = rep("b1", 2),
    pair_id = rep("pair1", 2),
    population = rep("PopA", 2),
    known_relationship = rep("parent_child", 2),
    loci_set = rep("set1", 2),
    tested_relationship = rep("parent_child", 2),
    tested_population = rep("PopA", 2),
    is_correct_pop = rep(TRUE, 2),
    combined_LR = c(100, 100),
    stringsAsFactors = FALSE
  )
  
  result <- calculate_ratio_stats(mock_correct_only)
  
  # No ratios to compute when there's no wrong population data
  expect_equal(nrow(result$combined_lrs_ratio), 0)
})

# =============================================================================
# TEST: Edge cases
# =============================================================================

test_that("calculate_summary_stats handles single observation", {
  mock_single <- data.frame(
    batch_id = "b1",
    pair_id = "p1",
    population = "PopA",
    known_relationship = "parent_child",
    loci_set = "set1",
    tested_relationship = "parent_child",
    tested_population = "PopA",
    is_correct_pop = TRUE,
    combined_LR = 100,
    stringsAsFactors = FALSE
  )
  
  result <- calculate_summary_stats(mock_single)
  
  expect_equal(result$n, 1)
  expect_equal(result$mean_LR, 100)
  expect_equal(result$median_LR, 100)
  expect_true(is.na(result$sd_LR))  # SD of single value is NA
})

test_that("calculate_cutoffs handles edge case with few unrelated pairs", {
  mock_few <- data.frame(
    batch_id = rep("b1", 3),
    pair_id = paste0("p", 1:3),
    population = rep("PopA", 3),
    known_relationship = rep("unrelated", 3),
    loci_set = rep("set1", 3),
    tested_relationship = rep("parent_child", 3),
    tested_population = rep("PopA", 3),
    is_correct_pop = rep(TRUE, 3),
    combined_LR = c(1, 2, 3),
    stringsAsFactors = FALSE
  )
  
  fp_rates <- c(1, 0.1, 0.01)
  result <- calculate_cutoffs(mock_few, fp_rates)
  
  expect_equal(result$n_unrelated, 3)
  # Cutoffs should still be computed
  expect_false(is.na(result$cutoff_1))
})
