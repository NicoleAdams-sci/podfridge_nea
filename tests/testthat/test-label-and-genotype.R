# Test label_and_genotype()
# This function maps alleles to letters and creates genotype pattern strings

library(testthat)

# Define the function locally for testing (copied from LR_kinship_utility_functions.R)
label_and_genotype <- function(alleles1, alleles2) {
  # Unique value sets
  shared_vals      <- intersect(alleles1, alleles2)
  unique_ind1_vals <- setdiff(alleles1, shared_vals)
  unique_ind2_vals <- setdiff(alleles2, c(shared_vals, unique_ind1_vals))
  
  # Deterministic letter assignment: shared first, then uniques of ind1,
  # then uniques of ind2.  Each category sorted so results never depend on
  # input order of alleles.
  allele_map <- list()
  next_letter <- 1
  
  for (v in sort(shared_vals)) {
    allele_map[[LETTERS[next_letter]]] <- v
    next_letter <- next_letter + 1
  }
  for (v in sort(unique_ind1_vals)) {
    allele_map[[LETTERS[next_letter]]] <- v
    next_letter <- next_letter + 1
  }
  for (v in sort(unique_ind2_vals)) {
    allele_map[[LETTERS[next_letter]]] <- v
    next_letter <- next_letter + 1
  }
  allele_map <- unlist(allele_map)      # named character vector
  
  # Translate each individual's alleles → letters, sort within genotype
  translate <- function(vec) {
    letters <- names(allele_map)[match(vec, allele_map)]
    paste(sort(letters), collapse = "")
  }
  geno1 <- translate(alleles1)          # "AA", "AB", etc.
  geno2 <- translate(alleles2)
  
  # Retain individual order, do NOT sort across genotypes
  genotype_match <- paste(geno1, geno2, sep = "-")
  
  list(
    allele_map     = allele_map,
    genotype_ind1  = geno1,
    genotype_ind2  = geno2,
    genotype_match = genotype_match
  )
}

# =============================================================================
# TEST: 2 shared alleles - identical genotypes
# =============================================================================

test_that("identical heterozygous genotypes produce AB-AB pattern", {
  result <- label_and_genotype(c("15", "16"), c("15", "16"))
  expect_equal(result$genotype_match, "AB-AB")
  expect_equal(result$genotype_ind1, "AB")
  expect_equal(result$genotype_ind2, "AB")
})

test_that("identical homozygous genotypes produce AA-AA pattern", {
  result <- label_and_genotype(c("15", "15"), c("15", "15"))
  expect_equal(result$genotype_match, "AA-AA")
  expect_equal(result$genotype_ind1, "AA")
  expect_equal(result$genotype_ind2, "AA")
})

# =============================================================================
# TEST: 1 shared allele patterns
# =============================================================================

test_that("homozygous vs heterozygous (sharing one) produces AA-AB pattern", {
  # ind1 is AA, ind2 is AB (shares A)
  result <- label_and_genotype(c("15", "15"), c("15", "16"))
  expect_equal(result$genotype_match, "AA-AB")
})

test_that("heterozygous vs homozygous (sharing one) produces AB-AA pattern", {
  # ind1 is AB, ind2 is AA (shares A)
  result <- label_and_genotype(c("15", "16"), c("15", "15"))
  expect_equal(result$genotype_match, "AB-AA")
})

test_that("heterozygous vs heterozygous (one shared) produces AB-AC pattern", {
  # ind1 is AB, ind2 is AC (shares A)
  result <- label_and_genotype(c("15", "16"), c("15", "17"))
  expect_equal(result$genotype_match, "AB-AC")
})

# =============================================================================
# TEST: 0 shared alleles
# =============================================================================

test_that("completely different genotypes produce AB-CD pattern", {
  result <- label_and_genotype(c("15", "16"), c("17", "18"))
  expect_equal(result$genotype_match, "AB-CD")
})

test_that("different homozygous genotypes produce AA-BB pattern", {
  result <- label_and_genotype(c("15", "15"), c("16", "16"))
  expect_equal(result$genotype_match, "AA-BB")
})

# =============================================================================
# TEST: Allele order independence
# =============================================================================

test_that("genotype pattern is independent of allele order within individual", {
  # (15, 16) and (16, 15) should give same genotype pattern
  result1 <- label_and_genotype(c("15", "16"), c("15", "17"))
  result2 <- label_and_genotype(c("16", "15"), c("17", "15"))
  expect_equal(result1$genotype_match, result2$genotype_match)
})

# =============================================================================
# TEST: Allele map correctness
# =============================================================================

test_that("allele map assigns shared alleles first", {
  result <- label_and_genotype(c("15", "16"), c("15", "17"))
  # "15" is shared, should be A
  # "16" is unique to ind1, should be B
  # "17" is unique to ind2, should be C
  expect_equal(result$allele_map[["A"]], "15")
  expect_equal(result$allele_map[["B"]], "16")
  expect_equal(result$allele_map[["C"]], "17")
})

test_that("allele map handles two shared alleles correctly", {
  result <- label_and_genotype(c("15", "16"), c("15", "16"))
  # Both shared, assigned in sorted order
  expect_equal(result$allele_map[["A"]], "15")
  expect_equal(result$allele_map[["B"]], "16")
})
