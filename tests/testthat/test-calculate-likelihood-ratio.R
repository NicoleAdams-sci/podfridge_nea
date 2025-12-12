# Test calculate_likelihood_ratio()
# This is the core LR formula function - pure math, no data dependencies

library(testthat)

# Define the function locally for testing (copied from LR_kinship_utility_functions.R)
calculate_likelihood_ratio <- function(shared_alleles,
                                       genotype_match,
                                       pA = NA_real_, pB = NA_real_,
                                       k0, k1, k2) {
  
  # Probability that the tested relationship shares 0 IBD alleles
  if (shared_alleles == 0) return(k0)
  
  if (shared_alleles == 1) {
    Rxp <- switch(genotype_match,
                  "AA-AB" = 2 * pA,
                  "AB-AA" = 2 * pA,
                  "AB-AC" = 4 * pA,
                  stop("Invalid genotype_match for 1 shared allele.")
    )
    return(k0 + (k1 / Rxp))
  }
  
  if (shared_alleles == 2) {
    if (genotype_match == "AA-AA") {
      Rxp <- pA
      Rxu <- pA^2
    } else if (genotype_match == "AB-AB") {
      Rxp <- (4 * pA * pB) / (pA + pB)
      Rxu <- 2 * pA * pB
    } else {
      stop("Invalid genotype_match for 2 shared alleles.")
    }
    return(k0 + (k1 / Rxp) + (k2 / Rxu))
  }
  
  stop("shared_alleles must be 0, 1, or 2.")
}

# =============================================================================
# Kinship coefficients for reference:
#   parent_child:   k0=0,    k1=1,    k2=0
#   full_siblings:  k0=0.25, k1=0.5,  k2=0.25
#   half_siblings:  k0=0.5,  k1=0.5,  k2=0
#   unrelated:      k0=1,    k1=0,    k2=0
# =============================================================================

# =============================================================================
# TEST: 0 shared alleles -> LR = k0
# =============================================================================

test_that("0 shared alleles returns k0 for parent-child", {
  # Parent-child: k0=0, k1=1, k2=0
  # If no shared alleles, LR = k0 = 0 (exclusion)
  lr <- calculate_likelihood_ratio(
    shared_alleles = 0,
    genotype_match = "AB-CD",  # doesn't matter for 0 shared
    pA = 0.25, pB = 0.30,      # doesn't matter for 0 shared
    k0 = 0, k1 = 1, k2 = 0
  )
  expect_equal(lr, 0)
})

test_that("0 shared alleles returns k0 for unrelated", {
  # Unrelated: k0=1, k1=0, k2=0
  lr <- calculate_likelihood_ratio(
    shared_alleles = 0,
    genotype_match = "AB-CD",
    pA = 0.25, pB = 0.30,
    k0 = 1, k1 = 0, k2 = 0
  )
  expect_equal(lr, 1)
})

test_that("0 shared alleles returns k0 for full siblings", {
  # Full siblings: k0=0.25
  lr <- calculate_likelihood_ratio(
    shared_alleles = 0,
    genotype_match = "AB-CD",
    pA = 0.25, pB = 0.30,
    k0 = 0.25, k1 = 0.5, k2 = 0.25
  )
  expect_equal(lr, 0.25)
})

# =============================================================================
# TEST: 1 shared allele - AA-AB pattern
# Formula: LR = k0 + k1/(2*pA)
# =============================================================================

test_that("AA-AB pattern calculates correctly for parent-child", {
  # Parent-child: k0=0, k1=1, k2=0
  # pA = 0.25
  # LR = 0 + 1/(2*0.25) = 1/0.5 = 2.0
  lr <- calculate_likelihood_ratio(
    shared_alleles = 1,
    genotype_match = "AA-AB",
    pA = 0.25, pB = NA,
    k0 = 0, k1 = 1, k2 = 0
  )
  expect_equal(lr, 2.0)
})

test_that("AA-AB pattern calculates correctly for full siblings", {
  # Full siblings: k0=0.25, k1=0.5, k2=0.25
  # pA = 0.25
  # LR = 0.25 + 0.5/(2*0.25) = 0.25 + 1.0 = 1.25
  lr <- calculate_likelihood_ratio(
    shared_alleles = 1,
    genotype_match = "AA-AB",
    pA = 0.25, pB = NA,
    k0 = 0.25, k1 = 0.5, k2 = 0.25
  )
  expect_equal(lr, 1.25)
})

# =============================================================================
# TEST: 1 shared allele - AB-AA pattern (symmetric to AA-AB)
# Formula: LR = k0 + k1/(2*pA)
# =============================================================================

test_that("AB-AA pattern calculates correctly for parent-child", {
  # Same formula as AA-AB
  lr <- calculate_likelihood_ratio(
    shared_alleles = 1,
    genotype_match = "AB-AA",
    pA = 0.25, pB = NA,
    k0 = 0, k1 = 1, k2 = 0
  )
  expect_equal(lr, 2.0)
})

# =============================================================================
# TEST: 1 shared allele - AB-AC pattern
# Formula: LR = k0 + k1/(4*pA)
# =============================================================================

test_that("AB-AC pattern calculates correctly for parent-child", {
  # Parent-child: k0=0, k1=1, k2=0
  # pA = 0.25
  # LR = 0 + 1/(4*0.25) = 1/1.0 = 1.0
  lr <- calculate_likelihood_ratio(
    shared_alleles = 1,
    genotype_match = "AB-AC",
    pA = 0.25, pB = NA,
    k0 = 0, k1 = 1, k2 = 0
  )
  expect_equal(lr, 1.0)
})

test_that("AB-AC pattern calculates correctly for half siblings", {
  # Half siblings: k0=0.5, k1=0.5, k2=0
  # pA = 0.20
  # LR = 0.5 + 0.5/(4*0.20) = 0.5 + 0.5/0.8 = 0.5 + 0.625 = 1.125
  lr <- calculate_likelihood_ratio(
    shared_alleles = 1,
    genotype_match = "AB-AC",
    pA = 0.20, pB = NA,
    k0 = 0.5, k1 = 0.5, k2 = 0
  )
  expect_equal(lr, 1.125)
})

# =============================================================================
# TEST: 2 shared alleles - AA-AA pattern
# Formula: LR = k0 + k1/pA + k2/pA^2
# =============================================================================

test_that("AA-AA pattern calculates correctly for parent-child", {
  # Parent-child: k0=0, k1=1, k2=0
  # pA = 0.25
  # LR = 0 + 1/0.25 + 0/0.0625 = 4.0
  lr <- calculate_likelihood_ratio(
    shared_alleles = 2,
    genotype_match = "AA-AA",
    pA = 0.25, pB = NA,
    k0 = 0, k1 = 1, k2 = 0
  )
  expect_equal(lr, 4.0)
})

test_that("AA-AA pattern calculates correctly for full siblings", {
  # Full siblings: k0=0.25, k1=0.5, k2=0.25
  # pA = 0.25
  # LR = 0.25 + 0.5/0.25 + 0.25/0.0625 = 0.25 + 2.0 + 4.0 = 6.25
  lr <- calculate_likelihood_ratio(
    shared_alleles = 2,
    genotype_match = "AA-AA",
    pA = 0.25, pB = NA,
    k0 = 0.25, k1 = 0.5, k2 = 0.25
  )
  expect_equal(lr, 6.25)
})

test_that("AA-AA pattern calculates correctly for unrelated", {
  # Unrelated: k0=1, k1=0, k2=0
  # LR = 1 + 0 + 0 = 1.0
  lr <- calculate_likelihood_ratio(
    shared_alleles = 2,
    genotype_match = "AA-AA",
    pA = 0.25, pB = NA,
    k0 = 1, k1 = 0, k2 = 0
  )
  expect_equal(lr, 1.0)
})

# =============================================================================
# TEST: 2 shared alleles - AB-AB pattern
# Formula: LR = k0 + k1*(pA+pB)/(4*pA*pB) + k2/(2*pA*pB)
# =============================================================================

test_that("AB-AB pattern calculates correctly for parent-child", {
  # Parent-child: k0=0, k1=1, k2=0
  # pA = 0.25, pB = 0.30
  # Rxp = 4*pA*pB/(pA+pB) = 4*0.25*0.30/(0.55) = 0.30/0.55 = 0.5454545...
  # LR = 0 + 1/Rxp + 0 = 1/0.5454545 = 1.8333...
  lr <- calculate_likelihood_ratio(
    shared_alleles = 2,
    genotype_match = "AB-AB",
    pA = 0.25, pB = 0.30,
    k0 = 0, k1 = 1, k2 = 0
  )
  expected <- 1 * (0.25 + 0.30) / (4 * 0.25 * 0.30)  # = 0.55 / 0.30 = 1.8333...
  expect_equal(lr, expected, tolerance = 1e-10)
})

test_that("AB-AB pattern calculates correctly for full siblings", {
  # Full siblings: k0=0.25, k1=0.5, k2=0.25
  # pA = 0.25, pB = 0.30
  # Rxp = 4*pA*pB/(pA+pB) = 0.30/0.55
  # Rxu = 2*pA*pB = 0.15
  # LR = 0.25 + 0.5/Rxp + 0.25/Rxu
  #    = 0.25 + 0.5*(0.55/0.30) + 0.25/0.15
  #    = 0.25 + 0.9166... + 1.6666...
  #    = 2.8333...
  lr <- calculate_likelihood_ratio(
    shared_alleles = 2,
    genotype_match = "AB-AB",
    pA = 0.25, pB = 0.30,
    k0 = 0.25, k1 = 0.5, k2 = 0.25
  )
  Rxp <- (4 * 0.25 * 0.30) / (0.25 + 0.30)
  Rxu <- 2 * 0.25 * 0.30
  expected <- 0.25 + 0.5/Rxp + 0.25/Rxu
  expect_equal(lr, expected, tolerance = 1e-10)
})

test_that("AB-AB pattern calculates correctly for unrelated", {
  # Unrelated: k0=1, k1=0, k2=0
  # LR = 1 + 0 + 0 = 1.0
  lr <- calculate_likelihood_ratio(
    shared_alleles = 2,
    genotype_match = "AB-AB",
    pA = 0.25, pB = 0.30,
    k0 = 1, k1 = 0, k2 = 0
  )
  expect_equal(lr, 1.0)
})

# =============================================================================
# TEST: Effect of allele frequency on LR
# =============================================================================

test_that("rarer alleles produce higher LRs for AA-AA parent-child", {
  # Rarer allele (pA = 0.05) should give higher LR than common allele (pA = 0.25)
  # For parent-child AA-AA: LR = k1/pA = 1/pA
  lr_rare <- calculate_likelihood_ratio(2, "AA-AA", pA = 0.05, pB = NA, k0 = 0, k1 = 1, k2 = 0)
  lr_common <- calculate_likelihood_ratio(2, "AA-AA", pA = 0.25, pB = NA, k0 = 0, k1 = 1, k2 = 0)
  
  expect_equal(lr_rare, 20.0)    # 1/0.05
  expect_equal(lr_common, 4.0)   # 1/0.25
  expect_gt(lr_rare, lr_common)
})

# =============================================================================
# TEST: Error handling
# =============================================================================

test_that("invalid genotype_match for 1 shared allele throws error", {
  expect_error(
    calculate_likelihood_ratio(1, "AB-AB", pA = 0.25, pB = 0.30, k0 = 0, k1 = 1, k2 = 0),
    "Invalid genotype_match"
  )
})

test_that("invalid genotype_match for 2 shared alleles throws error", {
  expect_error(
    calculate_likelihood_ratio(2, "AB-AC", pA = 0.25, pB = 0.30, k0 = 0, k1 = 1, k2 = 0),
    "Invalid genotype_match"
  )
})

test_that("invalid shared_alleles value throws error", {
  expect_error(
    calculate_likelihood_ratio(3, "AA-AA", pA = 0.25, pB = NA, k0 = 0, k1 = 1, k2 = 0),
    "shared_alleles must be 0, 1, or 2"
  )
})
