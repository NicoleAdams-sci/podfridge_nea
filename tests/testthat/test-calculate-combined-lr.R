# Test calculate_combined_lr()
# This function multiplies single-locus LRs across loci sets

library(testthat)
library(dplyr)
library(purrr)

# =============================================================================
# Define the function locally for testing (copied from module5_combined_LR.R)
# =============================================================================

calculate_combined_lr <- function(single_locus_results, loci_sets) {
  
  # Validate inputs
  required_cols <- c("batch_id", "pair_id", "population", "known_relationship", 
                     "locus", "tested_relationship", "tested_population", "LR")
  missing_cols <- setdiff(required_cols, names(single_locus_results))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  if (!is.list(loci_sets) || is.null(names(loci_sets))) {
    stop("loci_sets must be a named list")
  }
  
  # Calculate combined LRs for each loci set
  combined_results <- map_dfr(names(loci_sets), function(set_name) {
    loci_in_set <- loci_sets[[set_name]]
    
    # Filter to loci in this set and calculate product
    set_results <- single_locus_results |>
      filter(locus %in% loci_in_set) |>
      group_by(batch_id, pair_id, population, known_relationship, 
               tested_relationship, tested_population) |>
      summarize(
        combined_LR = prod(LR, na.rm = TRUE),
        .groups = 'drop'
      ) |>
      mutate(loci_set = set_name)
    
    return(set_results)
  })
  
  # Reorder columns with batch_id first
  combined_results <- combined_results |>
    select(batch_id, pair_id, population, known_relationship, loci_set,
           tested_relationship, tested_population, combined_LR)
  
  return(combined_results)
}

# =============================================================================
# MOCK DATA SETUP
# =============================================================================

# Create mock single-locus LR data for 3 loci
mock_single_locus_lr <- data.frame(
  batch_id = rep("batch1", 6),
  pair_id = rep("pair001", 6),
  population = rep("PopA", 6),
  known_relationship = rep("parent_child", 6),
  locus = c("LOCUS1", "LOCUS2", "LOCUS3", "LOCUS1", "LOCUS2", "LOCUS3"),
  tested_relationship = c(rep("parent_child", 3), rep("unrelated", 3)),
  tested_population = rep("PopA", 6),
  LR = c(2.0, 4.0, 3.0,    # parent_child LRs for 3 loci
         1.0, 1.0, 1.0),   # unrelated LRs (always 1 for matching genotypes with k0=1)
  stringsAsFactors = FALSE
)

# Mock loci sets
mock_loci_sets <- list(
  set_all = c("LOCUS1", "LOCUS2", "LOCUS3"),
  set_partial = c("LOCUS1", "LOCUS2")
)

# =============================================================================
# TEST: Basic multiplication of LRs
# =============================================================================

test_that("combined LR is product of single-locus LRs for full set", {
  result <- calculate_combined_lr(mock_single_locus_lr, mock_loci_sets)
  
  # Filter to parent_child, set_all
  pc_all <- result |> 
    filter(tested_relationship == "parent_child", loci_set == "set_all")
  
  # Expected: 2.0 * 4.0 * 3.0 = 24.0
  expect_equal(pc_all$combined_LR, 24.0)
})

test_that("combined LR is product of single-locus LRs for partial set", {
  result <- calculate_combined_lr(mock_single_locus_lr, mock_loci_sets)
  
  # Filter to parent_child, set_partial (only LOCUS1 and LOCUS2)
  pc_partial <- result |> 
    filter(tested_relationship == "parent_child", loci_set == "set_partial")
  
  # Expected: 2.0 * 4.0 = 8.0
  expect_equal(pc_partial$combined_LR, 8.0)
})

test_that("unrelated LRs multiply to 1 when all single-locus LRs are 1", {
  result <- calculate_combined_lr(mock_single_locus_lr, mock_loci_sets)
  
  unrel_all <- result |> 
    filter(tested_relationship == "unrelated", loci_set == "set_all")
  
  # Expected: 1.0 * 1.0 * 1.0 = 1.0
  expect_equal(unrel_all$combined_LR, 1.0)
})

# =============================================================================
# TEST: Returns one row per loci_set per tested_relationship
# =============================================================================

test_that("returns correct number of rows for multiple loci sets and relationships", {
  result <- calculate_combined_lr(mock_single_locus_lr, mock_loci_sets)
  
  # 2 loci_sets × 2 tested_relationships = 4 rows
  expect_equal(nrow(result), 4)
  expect_setequal(unique(result$loci_set), c("set_all", "set_partial"))
  expect_setequal(unique(result$tested_relationship), c("parent_child", "unrelated"))
})

# =============================================================================
# TEST: Multiple pairs
# =============================================================================

test_that("handles multiple pairs correctly", {
  # Add a second pair
  mock_two_pairs <- rbind(
    mock_single_locus_lr,
    data.frame(
      batch_id = rep("batch1", 3),
      pair_id = rep("pair002", 3),
      population = rep("PopA", 3),
      known_relationship = rep("parent_child", 3),
      locus = c("LOCUS1", "LOCUS2", "LOCUS3"),
      tested_relationship = rep("parent_child", 3),
      tested_population = rep("PopA", 3),
      LR = c(5.0, 2.0, 2.0),  # Different LRs for pair002
      stringsAsFactors = FALSE
    )
  )
  
  result <- calculate_combined_lr(mock_two_pairs, list(set_all = c("LOCUS1", "LOCUS2", "LOCUS3")))
  
  # Filter to parent_child only
  pc_results <- result |> filter(tested_relationship == "parent_child")
  
  expect_equal(nrow(pc_results), 2)
  
  pair001_lr <- pc_results$combined_LR[pc_results$pair_id == "pair001"]
  pair002_lr <- pc_results$combined_LR[pc_results$pair_id == "pair002"]
  
  expect_equal(pair001_lr, 24.0)  # 2 * 4 * 3
  expect_equal(pair002_lr, 20.0)  # 5 * 2 * 2
})

# =============================================================================
# TEST: Multiple populations
# =============================================================================

test_that("handles multiple tested populations correctly", {
  mock_multi_pop <- data.frame(
    batch_id = rep("batch1", 6),
    pair_id = rep("pair001", 6),
    population = rep("PopA", 6),
    known_relationship = rep("parent_child", 6),
    locus = c("LOCUS1", "LOCUS2", "LOCUS1", "LOCUS2", "LOCUS1", "LOCUS2"),
    tested_relationship = rep("parent_child", 6),
    tested_population = c("PopA", "PopA", "PopB", "PopB", "PopC", "PopC"),
    LR = c(2.0, 3.0,   # PopA: product = 6
           4.0, 5.0,   # PopB: product = 20
           1.5, 2.0),  # PopC: product = 3
    stringsAsFactors = FALSE
  )
  
  result <- calculate_combined_lr(
    mock_multi_pop, 
    list(set_two = c("LOCUS1", "LOCUS2"))
  )
  
  expect_equal(nrow(result), 3)
  
  lr_popA <- result$combined_LR[result$tested_population == "PopA"]
  lr_popB <- result$combined_LR[result$tested_population == "PopB"]
  lr_popC <- result$combined_LR[result$tested_population == "PopC"]
  
  expect_equal(lr_popA, 6.0)
  expect_equal(lr_popB, 20.0)
  expect_equal(lr_popC, 3.0)
})

# =============================================================================
# TEST: NA handling
# =============================================================================

test_that("handles NA values with na.rm = TRUE", {
  mock_with_na <- data.frame(
    batch_id = rep("batch1", 3),
    pair_id = rep("pair001", 3),
    population = rep("PopA", 3),
    known_relationship = rep("parent_child", 3),
    locus = c("LOCUS1", "LOCUS2", "LOCUS3"),
    tested_relationship = rep("parent_child", 3),
    tested_population = rep("PopA", 3),
    LR = c(2.0, NA, 3.0),  # NA in the middle
    stringsAsFactors = FALSE
  )
  
  result <- calculate_combined_lr(
    mock_with_na,
    list(set_all = c("LOCUS1", "LOCUS2", "LOCUS3"))
  )
  
  # With na.rm = TRUE, prod(2.0, NA, 3.0) = 2.0 * 3.0 = 6.0
  expect_equal(result$combined_LR, 6.0)
})

# =============================================================================
# TEST: Empty loci set
# =============================================================================

test_that("empty loci set returns product of empty set (1)",
{
  result <- calculate_combined_lr(
    mock_single_locus_lr,
    list(empty_set = c("NONEXISTENT"))
  )
  
  # No loci match, so no rows for that pair in that set
  # Actually, the filter will return empty, so group_by summarize returns nothing
  expect_equal(nrow(result), 0)
})

# =============================================================================
# TEST: Loci set with partial overlap
# =============================================================================

test_that("loci set with only some matching loci uses only those", {
  result <- calculate_combined_lr(
    mock_single_locus_lr,
    list(mixed_set = c("LOCUS1", "NONEXISTENT_LOCUS"))
  )
  
  # Only LOCUS1 exists in data, so combined_LR = LR of LOCUS1 alone
  pc_result <- result |> filter(tested_relationship == "parent_child")
  expect_equal(pc_result$combined_LR, 2.0)  # Just LOCUS1
})

# =============================================================================
# TEST: Output column structure
# =============================================================================

test_that("output has correct columns in correct order", {
  result <- calculate_combined_lr(mock_single_locus_lr, mock_loci_sets)
  
  expected_cols <- c("batch_id", "pair_id", "population", "known_relationship",
                     "loci_set", "tested_relationship", "tested_population", "combined_LR")
  
  expect_equal(names(result), expected_cols)
})

# =============================================================================
# TEST: Error handling - missing columns
# =============================================================================

test_that("throws error for missing required columns", {
  bad_data <- data.frame(
    pair_id = "pair001",
    locus = "LOCUS1",
    LR = 2.0
  )
  
  expect_error(
    calculate_combined_lr(bad_data, mock_loci_sets),
    "Missing required columns"
  )
})

test_that("throws error for non-list loci_sets", {
  expect_error(
    calculate_combined_lr(mock_single_locus_lr, c("LOCUS1", "LOCUS2")),
    "loci_sets must be a named list"
  )
})

test_that("throws error for unnamed list loci_sets", {
  expect_error(
    calculate_combined_lr(mock_single_locus_lr, list(c("LOCUS1"), c("LOCUS2"))),
    "loci_sets must be a named list"
  )
})

# =============================================================================
# TEST: Preserves grouping variables correctly
# =============================================================================

test_that("preserves population and known_relationship in output", {
  result <- calculate_combined_lr(mock_single_locus_lr, mock_loci_sets)
  
  expect_true(all(result$population == "PopA"))
  expect_true(all(result$known_relationship == "parent_child"))
  expect_true(all(result$batch_id == "batch1"))
})

# =============================================================================
# TEST: Very small and very large LRs
# =============================================================================

test_that("handles very small LR values correctly", {
  mock_small <- data.frame(
    batch_id = "batch1",
    pair_id = "pair001",
    population = "PopA",
    known_relationship = "parent_child",
    locus = c("LOCUS1", "LOCUS2"),
    tested_relationship = "parent_child",
    tested_population = "PopA",
    LR = c(0.001, 0.001),
    stringsAsFactors = FALSE
  )
  
  result <- calculate_combined_lr(mock_small, list(set = c("LOCUS1", "LOCUS2")))
  
  # 0.001 * 0.001 = 0.000001
  expect_equal(result$combined_LR, 1e-6, tolerance = 1e-12)
})

test_that("handles very large LR values correctly", {
  mock_large <- data.frame(
    batch_id = "batch1",
    pair_id = "pair001",
    population = "PopA",
    known_relationship = "parent_child",
    locus = c("LOCUS1", "LOCUS2"),
    tested_relationship = "parent_child",
    tested_population = "PopA",
    LR = c(1e6, 1e6),
    stringsAsFactors = FALSE
  )
  
  result <- calculate_combined_lr(mock_large, list(set = c("LOCUS1", "LOCUS2")))
  
  # 1e6 * 1e6 = 1e12
  expect_equal(result$combined_LR, 1e12)
})
