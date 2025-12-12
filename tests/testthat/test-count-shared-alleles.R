# Test count_shared_alleles()
# This is the smallest unit - counts how many alleles are shared between two individuals

library(testthat)

# Define the function locally for testing (copied from LR_kinship_utility_functions.R)
# This avoids loading the entire module which tries to read data files
count_shared_alleles <- function(alleles1, alleles2) {
  alleles1 <- as.vector(as.character(alleles1))
  alleles2 <- as.vector(as.character(alleles2))
  
  alleles1 <- alleles1[!is.na(alleles1)]
  alleles2 <- alleles2[!is.na(alleles2)]
  
  if (length(alleles1) == 0 || length(alleles2) == 0) {
    return(0)
  }
  
  tbl1 <- table(alleles1)
  tbl2 <- table(alleles2)
  common <- intersect(names(tbl1), names(tbl2))
  
  if (length(common) == 0) {
    return(0)
  }
  
  sum(mapply(function(a) min(tbl1[[a]], tbl2[[a]]), common))
}

# =============================================================================
# TEST: 2 shared alleles
# =============================================================================

test_that("identical heterozygous genotypes share 2 alleles", {
  # AB vs AB -> 2 shared

  expect_equal(count_shared_alleles(c("15", "16"), c("15", "16")), 2)
  expect_equal(count_shared_alleles(c("12", "14"), c("14", "12")), 2)  # order shouldn't matter
})

test_that("identical homozygous genotypes share 2 alleles", {
  # AA vs AA -> 2 shared
  expect_equal(count_shared_alleles(c("15", "15"), c("15", "15")), 2)
  expect_equal(count_shared_alleles(c("9", "9"), c("9", "9")), 2)
})

# =============================================================================
# TEST: 1 shared allele
# =============================================================================

test_that("heterozygous vs heterozygous with one common allele shares 1", {

  # AB vs AC -> 1 shared (A)
  expect_equal(count_shared_alleles(c("15", "16"), c("15", "17")), 1)
  expect_equal(count_shared_alleles(c("12", "14"), c("13", "14")), 1)
})
  
test_that("homozygous vs heterozygous with common allele shares 1", {
  # AA vs AB -> 1 shared
  expect_equal(count_shared_alleles(c("15", "15"), c("15", "16")), 1)
  # AB vs AA -> 1 shared  
  expect_equal(count_shared_alleles(c("15", "16"), c("15", "15")), 1)
})

# =============================================================================
# TEST: 0 shared alleles
# =============================================================================

test_that("completely different genotypes share 0 alleles", {
  # AB vs CD -> 0 shared
  expect_equal(count_shared_alleles(c("15", "16"), c("17", "18")), 0)
  expect_equal(count_shared_alleles(c("9", "10"), c("11", "12")), 0)
})

test_that("homozygous vs different homozygous shares 0 alleles", {
  # AA vs BB -> 0 shared
  expect_equal(count_shared_alleles(c("15", "15"), c("16", "16")), 0)
})

# =============================================================================
# TEST: Edge cases
# =============================================================================

test_that("handles numeric input by converting to character", {
  expect_equal(count_shared_alleles(c(15, 16), c(15, 16)), 2)
  expect_equal(count_shared_alleles(c(15, 16), c(17, 18)), 0)
})

test_that("handles NA values gracefully", {
  expect_equal(count_shared_alleles(c(NA, NA), c("15", "16")), 0)
  expect_equal(count_shared_alleles(c("15", NA), c("15", "16")), 1)
})

test_that("handles empty vectors", {
  expect_equal(count_shared_alleles(c(), c("15", "16")), 0)
  expect_equal(count_shared_alleles(c("15", "16"), c()), 0)
})
