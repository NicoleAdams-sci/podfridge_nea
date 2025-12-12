# Test LR calculations for ALL relationship types across ALL genotype patterns
# This ensures comprehensive coverage of the likelihood ratio formulas
#
# Relationships tested:
# - parent_child:   k0=0,    k1=1,     k2=0
# - full_siblings:  k0=0.25, k1=0.5,   k2=0.25
# - half_siblings:  k0=0.5,  k1=0.5,   k2=0     (also covers avuncular, grandparent)
# - cousins:        k0=0.875,k1=0.125, k2=0
# - unrelated:      k0=1,    k1=0,     k2=0
#
# Genotype patterns tested:
# - 0 shared: AB-CD (any pattern)
# - 1 shared: AA-AB, AB-AA, AB-AC
# - 2 shared: AA-AA, AB-AB

library(testthat)

# =============================================================================
# Define the function locally (copied from LR_kinship_utility_functions.R)
# =============================================================================

calculate_likelihood_ratio <- function(shared_alleles,
                                       genotype_match,
                                       pA = NA_real_, pB = NA_real_,
                                       k0, k1, k2) {
  
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
# Define kinship coefficients for all relationships
# =============================================================================

KINSHIP_COEFFS <- list(
  parent_child  = list(k0 = 0,     k1 = 1,     k2 = 0),
  full_siblings = list(k0 = 0.25,  k1 = 0.5,   k2 = 0.25),
  half_siblings = list(k0 = 0.5,   k1 = 0.5,   k2 = 0),
  cousins       = list(k0 = 0.875, k1 = 0.125, k2 = 0),
  unrelated     = list(k0 = 1,     k1 = 0,     k2 = 0)
)

# Test allele frequencies
pA <- 0.1   # 10% frequency
pB <- 0.2   # 20% frequency

# =============================================================================
# PARENT-CHILD (k0=0, k1=1, k2=0)
# =============================================================================

test_that("parent_child: 0 shared alleles gives LR = 0", {
  k <- KINSHIP_COEFFS$parent_child
  # LR = k0 = 0
  expect_equal(
    calculate_likelihood_ratio(0, "AB-CD", pA, pB, k$k0, k$k1, k$k2),
    0
  )
})

test_that("parent_child: 1 shared AA-AB gives LR = k0 + k1/(2*pA)", {
  k <- KINSHIP_COEFFS$parent_child
  # LR = 0 + 1/(2*0.1) = 1/0.2 = 5
  expected <- 0 + 1 / (2 * pA)  # = 5
  expect_equal(
    calculate_likelihood_ratio(1, "AA-AB", pA, pB, k$k0, k$k1, k$k2),
    expected
  )
})

test_that("parent_child: 1 shared AB-AA gives LR = k0 + k1/(2*pA)", {
  k <- KINSHIP_COEFFS$parent_child
  # Same formula as AA-AB
  expected <- 0 + 1 / (2 * pA)  # = 5
  expect_equal(
    calculate_likelihood_ratio(1, "AB-AA", pA, pB, k$k0, k$k1, k$k2),
    expected
  )
})

test_that("parent_child: 1 shared AB-AC gives LR = k0 + k1/(4*pA)", {
  k <- KINSHIP_COEFFS$parent_child
  # LR = 0 + 1/(4*0.1) = 1/0.4 = 2.5
  expected <- 0 + 1 / (4 * pA)  # = 2.5
  expect_equal(
    calculate_likelihood_ratio(1, "AB-AC", pA, pB, k$k0, k$k1, k$k2),
    expected
  )
})

test_that("parent_child: 2 shared AA-AA gives LR = k0 + k1/pA + k2/pA^2", {

  k <- KINSHIP_COEFFS$parent_child
  # LR = 0 + 1/0.1 + 0/0.01 = 10
  expected <- 0 + 1 / pA + 0 / (pA^2)  # = 10
  expect_equal(
    calculate_likelihood_ratio(2, "AA-AA", pA, pB, k$k0, k$k1, k$k2),
    expected
  )
})

test_that("parent_child: 2 shared AB-AB gives correct LR", {
  k <- KINSHIP_COEFFS$parent_child
  # Rxp = 4*pA*pB/(pA+pB) = 4*0.1*0.2/0.3 = 0.08/0.3 = 0.2667
  # Rxu = 2*pA*pB = 2*0.1*0.2 = 0.04
  # LR = 0 + 1/0.2667 + 0/0.04 = 3.75
  Rxp <- (4 * pA * pB) / (pA + pB)
  expected <- 0 + 1 / Rxp + 0  # = 3.75
  expect_equal(
    calculate_likelihood_ratio(2, "AB-AB", pA, pB, k$k0, k$k1, k$k2),
    expected,
    tolerance = 1e-10
  )
})

# =============================================================================
# FULL SIBLINGS (k0=0.25, k1=0.5, k2=0.25)
# =============================================================================

test_that("full_siblings: 0 shared alleles gives LR = 0.25", {
  k <- KINSHIP_COEFFS$full_siblings
  # LR = k0 = 0.25
  expect_equal(
    calculate_likelihood_ratio(0, "AB-CD", pA, pB, k$k0, k$k1, k$k2),
    0.25
  )
})

test_that("full_siblings: 1 shared AA-AB gives LR = k0 + k1/(2*pA)", {
  k <- KINSHIP_COEFFS$full_siblings
  # LR = 0.25 + 0.5/(2*0.1) = 0.25 + 2.5 = 2.75
  expected <- 0.25 + 0.5 / (2 * pA)  # = 2.75
  expect_equal(
    calculate_likelihood_ratio(1, "AA-AB", pA, pB, k$k0, k$k1, k$k2),
    expected
  )
})

test_that("full_siblings: 1 shared AB-AA gives LR = k0 + k1/(2*pA)", {
  k <- KINSHIP_COEFFS$full_siblings
  expected <- 0.25 + 0.5 / (2 * pA)  # = 2.75
  expect_equal(
    calculate_likelihood_ratio(1, "AB-AA", pA, pB, k$k0, k$k1, k$k2),
    expected
  )
})

test_that("full_siblings: 1 shared AB-AC gives LR = k0 + k1/(4*pA)", {
  k <- KINSHIP_COEFFS$full_siblings
  # LR = 0.25 + 0.5/(4*0.1) = 0.25 + 1.25 = 1.5
  expected <- 0.25 + 0.5 / (4 * pA)  # = 1.5
  expect_equal(
    calculate_likelihood_ratio(1, "AB-AC", pA, pB, k$k0, k$k1, k$k2),
    expected
  )
})

test_that("full_siblings: 2 shared AA-AA gives LR = k0 + k1/pA + k2/pA^2", {
  k <- KINSHIP_COEFFS$full_siblings
  # LR = 0.25 + 0.5/0.1 + 0.25/0.01 = 0.25 + 5 + 25 = 30.25
  expected <- 0.25 + 0.5 / pA + 0.25 / (pA^2)  # = 30.25
  expect_equal(
    calculate_likelihood_ratio(2, "AA-AA", pA, pB, k$k0, k$k1, k$k2),
    expected
  )
})

test_that("full_siblings: 2 shared AB-AB gives correct LR", {
  k <- KINSHIP_COEFFS$full_siblings
  # Rxp = 4*pA*pB/(pA+pB) = 0.08/0.3 = 0.2667
  # Rxu = 2*pA*pB = 0.04
  # LR = 0.25 + 0.5/0.2667 + 0.25/0.04 = 0.25 + 1.875 + 6.25 = 8.375
  Rxp <- (4 * pA * pB) / (pA + pB)
  Rxu <- 2 * pA * pB
  expected <- 0.25 + 0.5 / Rxp + 0.25 / Rxu
  expect_equal(
    calculate_likelihood_ratio(2, "AB-AB", pA, pB, k$k0, k$k1, k$k2),
    expected,
    tolerance = 1e-10
  )
})

# =============================================================================
# HALF SIBLINGS (k0=0.5, k1=0.5, k2=0) - also covers avuncular, grandparent
# =============================================================================

test_that("half_siblings: 0 shared alleles gives LR = 0.5", {
  k <- KINSHIP_COEFFS$half_siblings
  # LR = k0 = 0.5
  expect_equal(
    calculate_likelihood_ratio(0, "AB-CD", pA, pB, k$k0, k$k1, k$k2),
    0.5
  )
})

test_that("half_siblings: 1 shared AA-AB gives LR = k0 + k1/(2*pA)", {
  k <- KINSHIP_COEFFS$half_siblings
  # LR = 0.5 + 0.5/(2*0.1) = 0.5 + 2.5 = 3
  expected <- 0.5 + 0.5 / (2 * pA)  # = 3
  expect_equal(
    calculate_likelihood_ratio(1, "AA-AB", pA, pB, k$k0, k$k1, k$k2),
    expected
  )
})

test_that("half_siblings: 1 shared AB-AA gives LR = k0 + k1/(2*pA)", {
  k <- KINSHIP_COEFFS$half_siblings
  expected <- 0.5 + 0.5 / (2 * pA)  # = 3
  expect_equal(
    calculate_likelihood_ratio(1, "AB-AA", pA, pB, k$k0, k$k1, k$k2),
    expected
  )
})

test_that("half_siblings: 1 shared AB-AC gives LR = k0 + k1/(4*pA)", {
  k <- KINSHIP_COEFFS$half_siblings
  # LR = 0.5 + 0.5/(4*0.1) = 0.5 + 1.25 = 1.75
  expected <- 0.5 + 0.5 / (4 * pA)  # = 1.75
  expect_equal(
    calculate_likelihood_ratio(1, "AB-AC", pA, pB, k$k0, k$k1, k$k2),
    expected
  )
})

test_that("half_siblings: 2 shared AA-AA gives LR = k0 + k1/pA (k2=0)", {
  k <- KINSHIP_COEFFS$half_siblings
  # LR = 0.5 + 0.5/0.1 + 0/0.01 = 0.5 + 5 = 5.5
  # Note: k2=0 so no k2 term!
  expected <- 0.5 + 0.5 / pA + 0  # = 5.5
  expect_equal(
    calculate_likelihood_ratio(2, "AA-AA", pA, pB, k$k0, k$k1, k$k2),
    expected
  )
})

test_that("half_siblings: 2 shared AB-AB gives correct LR (k2=0)", {
  k <- KINSHIP_COEFFS$half_siblings
  # Rxp = 4*pA*pB/(pA+pB) = 0.08/0.3 = 0.2667
  # Rxu = 2*pA*pB = 0.04
  # LR = 0.5 + 0.5/0.2667 + 0/0.04 = 0.5 + 1.875 = 2.375
  # Note: k2=0 so no k2 term!
  Rxp <- (4 * pA * pB) / (pA + pB)
  expected <- 0.5 + 0.5 / Rxp + 0
  expect_equal(
    calculate_likelihood_ratio(2, "AB-AB", pA, pB, k$k0, k$k1, k$k2),
    expected,
    tolerance = 1e-10
  )
})

# =============================================================================
# COUSINS (k0=0.875, k1=0.125, k2=0)
# =============================================================================

test_that("cousins: 0 shared alleles gives LR = 0.875", {
  k <- KINSHIP_COEFFS$cousins
  # LR = k0 = 0.875
  expect_equal(
    calculate_likelihood_ratio(0, "AB-CD", pA, pB, k$k0, k$k1, k$k2),
    0.875
  )
})

test_that("cousins: 1 shared AA-AB gives LR = k0 + k1/(2*pA)", {
  k <- KINSHIP_COEFFS$cousins
  # LR = 0.875 + 0.125/(2*0.1) = 0.875 + 0.625 = 1.5
  expected <- 0.875 + 0.125 / (2 * pA)  # = 1.5
  expect_equal(
    calculate_likelihood_ratio(1, "AA-AB", pA, pB, k$k0, k$k1, k$k2),
    expected
  )
})

test_that("cousins: 1 shared AB-AA gives LR = k0 + k1/(2*pA)", {
  k <- KINSHIP_COEFFS$cousins
  expected <- 0.875 + 0.125 / (2 * pA)  # = 1.5
  expect_equal(
    calculate_likelihood_ratio(1, "AB-AA", pA, pB, k$k0, k$k1, k$k2),
    expected
  )
})

test_that("cousins: 1 shared AB-AC gives LR = k0 + k1/(4*pA)", {
  k <- KINSHIP_COEFFS$cousins
  # LR = 0.875 + 0.125/(4*0.1) = 0.875 + 0.3125 = 1.1875
  expected <- 0.875 + 0.125 / (4 * pA)  # = 1.1875
  expect_equal(
    calculate_likelihood_ratio(1, "AB-AC", pA, pB, k$k0, k$k1, k$k2),
    expected
  )
})

test_that("cousins: 2 shared AA-AA gives LR = k0 + k1/pA (k2=0)", {
  k <- KINSHIP_COEFFS$cousins
  # LR = 0.875 + 0.125/0.1 + 0 = 0.875 + 1.25 = 2.125
  expected <- 0.875 + 0.125 / pA + 0  # = 2.125
  expect_equal(
    calculate_likelihood_ratio(2, "AA-AA", pA, pB, k$k0, k$k1, k$k2),
    expected
  )
})

test_that("cousins: 2 shared AB-AB gives correct LR (k2=0)", {
  k <- KINSHIP_COEFFS$cousins
  # Rxp = 4*pA*pB/(pA+pB) = 0.08/0.3 = 0.2667
  # LR = 0.875 + 0.125/0.2667 + 0 = 0.875 + 0.46875 = 1.34375
  Rxp <- (4 * pA * pB) / (pA + pB)
  expected <- 0.875 + 0.125 / Rxp + 0
  expect_equal(
    calculate_likelihood_ratio(2, "AB-AB", pA, pB, k$k0, k$k1, k$k2),
    expected,
    tolerance = 1e-10
  )
})

# =============================================================================
# UNRELATED (k0=1, k1=0, k2=0)
# =============================================================================

test_that("unrelated: 0 shared alleles gives LR = 1", {
  k <- KINSHIP_COEFFS$unrelated
  # LR = k0 = 1
  expect_equal(
    calculate_likelihood_ratio(0, "AB-CD", pA, pB, k$k0, k$k1, k$k2),
    1
  )
})

test_that("unrelated: 1 shared AA-AB gives LR = 1 (k1=0)", {
  k <- KINSHIP_COEFFS$unrelated
  # LR = 1 + 0/(2*0.1) = 1
  expected <- 1 + 0 / (2 * pA)  # = 1
  expect_equal(
    calculate_likelihood_ratio(1, "AA-AB", pA, pB, k$k0, k$k1, k$k2),
    expected
  )
})

test_that("unrelated: 1 shared AB-AA gives LR = 1 (k1=0)", {
  k <- KINSHIP_COEFFS$unrelated
  expected <- 1 + 0 / (2 * pA)  # = 1
  expect_equal(
    calculate_likelihood_ratio(1, "AB-AA", pA, pB, k$k0, k$k1, k$k2),
    expected
  )
})

test_that("unrelated: 1 shared AB-AC gives LR = 1 (k1=0)", {
  k <- KINSHIP_COEFFS$unrelated
  # LR = 1 + 0/(4*0.1) = 1
  expected <- 1 + 0 / (4 * pA)  # = 1
  expect_equal(
    calculate_likelihood_ratio(1, "AB-AC", pA, pB, k$k0, k$k1, k$k2),
    expected
  )
})

test_that("unrelated: 2 shared AA-AA gives LR = 1 (k1=0, k2=0)", {
  k <- KINSHIP_COEFFS$unrelated
  # LR = 1 + 0/0.1 + 0/0.01 = 1
  expected <- 1  # all k1 and k2 terms are 0
  expect_equal(
    calculate_likelihood_ratio(2, "AA-AA", pA, pB, k$k0, k$k1, k$k2),
    expected
  )
})

test_that("unrelated: 2 shared AB-AB gives LR = 1 (k1=0, k2=0)", {
  k <- KINSHIP_COEFFS$unrelated
  # LR = 1 + 0 + 0 = 1
  expected <- 1
  expect_equal(
    calculate_likelihood_ratio(2, "AB-AB", pA, pB, k$k0, k$k1, k$k2),
    expected
  )
})

# =============================================================================
# CROSS-RELATIONSHIP COMPARISONS
# Verify that LRs follow expected patterns across relationships
# =============================================================================

test_that("parent_child has higher LR than full_siblings for 1 shared (same pA)", {
  pA_test <- 0.1
  
  pc_lr <- calculate_likelihood_ratio(1, "AA-AB", pA_test, NA, 0, 1, 0)
  fs_lr <- calculate_likelihood_ratio(1, "AA-AB", pA_test, NA, 0.25, 0.5, 0.25)
  
  # PC: 0 + 1/0.2 = 5
  # FS: 0.25 + 0.5/0.2 = 2.75
  expect_gt(pc_lr, fs_lr)
})

test_that("full_siblings has higher LR than half_siblings for 2 shared AA-AA", {
  pA_test <- 0.1
  
  # FS has k2=0.25, HS has k2=0, so FS gets a boost from k2 term
  fs_lr <- calculate_likelihood_ratio(2, "AA-AA", pA_test, NA, 0.25, 0.5, 0.25)
  hs_lr <- calculate_likelihood_ratio(2, "AA-AA", pA_test, NA, 0.5, 0.5, 0)
  
  # FS: 0.25 + 5 + 25 = 30.25
  # HS: 0.5 + 5 + 0 = 5.5
  expect_gt(fs_lr, hs_lr)
})

test_that("half_siblings has higher LR than cousins for same genotype", {
  pA_test <- 0.1
  
  hs_lr <- calculate_likelihood_ratio(1, "AA-AB", pA_test, NA, 0.5, 0.5, 0)
  co_lr <- calculate_likelihood_ratio(1, "AA-AB", pA_test, NA, 0.875, 0.125, 0)
  
  # HS: 0.5 + 2.5 = 3
  # CO: 0.875 + 0.625 = 1.5
  expect_gt(hs_lr, co_lr)
})

test_that("unrelated always gives LR = 1 regardless of shared alleles", {
  k <- KINSHIP_COEFFS$unrelated
  
  lr_0 <- calculate_likelihood_ratio(0, "AB-CD", 0.1, 0.2, k$k0, k$k1, k$k2)
  lr_1 <- calculate_likelihood_ratio(1, "AA-AB", 0.1, NA, k$k0, k$k1, k$k2)
  lr_2 <- calculate_likelihood_ratio(2, "AB-AB", 0.1, 0.2, k$k0, k$k1, k$k2)
  
  expect_equal(lr_0, 1)
  expect_equal(lr_1, 1)
  expect_equal(lr_2, 1)
})

# =============================================================================
# ALLELE FREQUENCY EFFECTS
# Verify that rarer alleles produce higher LRs (more discriminating)
# =============================================================================

test_that("rarer shared allele produces higher LR for parent_child", {
  k <- KINSHIP_COEFFS$parent_child
  
  common_lr <- calculate_likelihood_ratio(1, "AA-AB", 0.3, NA, k$k0, k$k1, k$k2)
  rare_lr <- calculate_likelihood_ratio(1, "AA-AB", 0.05, NA, k$k0, k$k1, k$k2)
  
  # Common (pA=0.3): 0 + 1/0.6 = 1.67
  # Rare (pA=0.05): 0 + 1/0.1 = 10
  expect_gt(rare_lr, common_lr)
})

test_that("rarer shared allele produces higher LR for half_siblings", {
  k <- KINSHIP_COEFFS$half_siblings
  
  common_lr <- calculate_likelihood_ratio(1, "AA-AB", 0.3, NA, k$k0, k$k1, k$k2)
  rare_lr <- calculate_likelihood_ratio(1, "AA-AB", 0.05, NA, k$k0, k$k1, k$k2)
  
  # Common (pA=0.3): 0.5 + 0.5/0.6 = 1.33
  # Rare (pA=0.05): 0.5 + 0.5/0.1 = 5.5
  expect_gt(rare_lr, common_lr)
})

test_that("k2 term amplifies LR for rare alleles in full_siblings", {
  k <- KINSHIP_COEFFS$full_siblings
  
  # For AA-AA pattern, k2 term is k2/pA^2, which grows rapidly as pA decreases
  common_lr <- calculate_likelihood_ratio(2, "AA-AA", 0.3, NA, k$k0, k$k1, k$k2)
  rare_lr <- calculate_likelihood_ratio(2, "AA-AA", 0.05, NA, k$k0, k$k1, k$k2)
  
  # Common: 0.25 + 0.5/0.3 + 0.25/0.09 = 0.25 + 1.67 + 2.78 = 4.7
  # Rare: 0.25 + 0.5/0.05 + 0.25/0.0025 = 0.25 + 10 + 100 = 110.25
  expect_gt(rare_lr, common_lr)
  expect_gt(rare_lr / common_lr, 10)  # Much larger difference due to k2 term
})
