# Integration Tests for the Kinship Calculation Pipeline
# 
# These tests verify that the full kinship_calculation() function works correctly
# end-to-end, using mock data that mirrors the real data structure.
#
# Unlike unit tests that test functions in isolation, integration tests verify:
# 1. Kinship coefficient lookup from the kinship_matrix
# 2. Allele frequency lookup from allele_frequency_data
# 3. Correct LR calculation combining all components
# 4. Proper handling of multiple tested populations
# 5. Edge cases like missing alleles or frequencies

library(testthat)
library(dplyr)
library(purrr)
library(tibble)

# =============================================================================
# Define ALL functions from the pipeline (copied from LR_kinship_utility_functions.R)
# =============================================================================

FALLBACK_FREQ <- 0.001  # Same as in the source file

count_shared_alleles <- function(alleles1, alleles2) {
  alleles1 <- as.vector(as.character(alleles1))
  alleles2 <- as.vector(as.character(alleles2))
  
  alleles1 <- alleles1[!is.na(alleles1)]
  alleles2 <- alleles2[!is.na(alleles2)]
  
  if (length(alleles1) == 0 || length(alleles2) == 0) return(0)
  
  tbl1 <- table(alleles1)
  tbl2 <- table(alleles2)
  common <- intersect(names(tbl1), names(tbl2))
  
  if (length(common) == 0) return(0)
  sum(mapply(function(a) min(tbl1[[a]], tbl2[[a]]), common))
}

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

kinship_calculation <- function(row,
                                allele_frequency_data,
                                kinship_matrix,
                                tested_relationship,
                                tested_populations) {
  
  if (anyNA(row[c("ind1_allele1", "ind1_allele2",
                  "ind2_allele1", "ind2_allele2")])) {
    out <- tibble::tibble(
      !!!row,
      shared_alleles      = NA_integer_,
      genotype_match      = NA_character_,
      tested_relationship = tested_relationship,
      tested_population   = tested_populations,
      LR                  = NA_real_
    )
    return(out)
  }
  
  alleles_ind1 <- as.character(c(row$ind1_allele1, row$ind1_allele2))
  alleles_ind2 <- as.character(c(row$ind2_allele1, row$ind2_allele2))
  
  shared_alleles <- count_shared_alleles(alleles_ind1, alleles_ind2)
  lab <- label_and_genotype(alleles_ind1, alleles_ind2)
  
  kin_vals <- kinship_matrix |>
    filter(relationship_type == tested_relationship)
  
  if (nrow(kin_vals) != 1L)
    stop("tested_relationship not found or duplicated in kinship_matrix")
  
  k0 <- kin_vals$k0
  k1 <- kin_vals$k1
  k2 <- kin_vals$k2
  
  fetch_freq <- function(freq_df, target_allele) {
    if (is.na(target_allele)) return(NA_real_)
    idx <- match(target_allele, freq_df$allele)
    if (is.na(idx)) return(FALLBACK_FREQ)
    f <- freq_df$frequency[[idx]]
    if (is.na(f) || f == 0) FALLBACK_FREQ else f
  }
  
  rows_out <- map_dfr(tested_populations, function(pop) {
    
    pop_freqs <- allele_frequency_data |>
      filter(population == pop, marker == row$locus)
    
    if (nrow(pop_freqs) == 0L) {
      LR_val <- NA_real_
    } else {
      A_val <- if ("A" %in% names(lab$allele_map)) lab$allele_map[["A"]] else NA_character_
      B_val <- if ("B" %in% names(lab$allele_map)) lab$allele_map[["B"]] else NA_character_
      
      pA <- fetch_freq(pop_freqs, A_val)
      pB <- fetch_freq(pop_freqs, B_val)
      
      if (is.na(pA)) {
        LR_val <- NA_real_
      } else if (shared_alleles == 2 &&
                 lab$genotype_match == "AB-AB" &&
                 is.na(pB)) {
        LR_val <- NA_real_
      } else {
        LR_val <- calculate_likelihood_ratio(
          shared_alleles,
          lab$genotype_match,
          pA, pB,
          k0, k1, k2
        )
      }
    }
    
    tibble::tibble(
      !!!row,
      shared_alleles      = shared_alleles,
      genotype_match      = lab$genotype_match,
      tested_relationship = tested_relationship,
      tested_population   = pop,
      LR                  = LR_val
    )
  })
  
  rows_out
}

# =============================================================================
# MOCK DATA SETUP - Mirrors real data structure
# =============================================================================

# Kinship coefficients (matches data/kinship_coefficients.csv)
mock_kinship_matrix <- data.frame(
  relationship_type = c("parent_child", "full_siblings", "half_siblings", "cousins", "unrelated"),
  k0 = c(0, 0.25, 0.5, 0.875, 1),
  k1 = c(1, 0.5, 0.5, 0.125, 0),
  k2 = c(0, 0.25, 0, 0, 0),
  stringsAsFactors = FALSE
)

# Allele frequencies for two populations
mock_allele_freq <- data.frame(
  population = c(rep("PopA", 6), rep("PopB", 6)),
  marker = rep(c("D3S1358", "D3S1358", "D3S1358", "vWA", "vWA", "vWA"), 2),
  allele = rep(c("15", "16", "17", "14", "15", "16"), 2),
  frequency = c(
    # PopA frequencies
    0.10, 0.20, 0.30,  # D3S1358: allele 15=10%, 16=20%, 17=30%
    0.15, 0.25, 0.35,  # vWA: allele 14=15%, 15=25%, 16=35%
    # PopB frequencies (different from PopA)
    0.25, 0.15, 0.20,  # D3S1358: allele 15=25%, 16=15%, 17=20%
    0.20, 0.30, 0.25   # vWA: allele 14=20%, 15=30%, 16=25%
  ),
  stringsAsFactors = FALSE
)

# Create a genotype row (mimics a row from paired genotype data)
create_genotype_row <- function(locus, ind1_a1, ind1_a2, ind2_a1, ind2_a2,
                                 pair_id = "pair1", batch_id = "batch1",
                                 population = "PopA", known_rel = "parent_child") {
  tibble::tibble(
    batch_id = batch_id,
    pair_id = pair_id,
    population = population,
    known_relationship = known_rel,
    locus = locus,
    ind1_allele1 = ind1_a1,
    ind1_allele2 = ind1_a2,
    ind2_allele1 = ind2_a1,
    ind2_allele2 = ind2_a2
  )
}

# =============================================================================
# TEST: Kinship coefficient lookup
# =============================================================================

test_that("kinship_calculation correctly looks up parent_child coefficients", {
  row <- create_genotype_row("D3S1358", "15", "16", "15", "17")
  
  result <- kinship_calculation(
    row,
    mock_allele_freq,
    mock_kinship_matrix,
    tested_relationship = "parent_child",
    tested_populations = "PopA"
  )
  
  # For AB-AC pattern with 1 shared allele:
  # LR = k0 + k1/(4*pA) = 0 + 1/(4*0.1) = 2.5
  expected_lr <- 0 + 1 / (4 * 0.10)  # pA for allele "15" in PopA
  
  expect_equal(result$LR, expected_lr)
  expect_equal(result$shared_alleles, 1)
  expect_equal(result$genotype_match, "AB-AC")
})

test_that("kinship_calculation correctly looks up half_siblings coefficients", {
  row <- create_genotype_row("D3S1358", "15", "16", "15", "17")
  
  result <- kinship_calculation(
    row,
    mock_allele_freq,
    mock_kinship_matrix,
    tested_relationship = "half_siblings",
    tested_populations = "PopA"
  )
  
  # For AB-AC pattern with 1 shared allele:
  # LR = k0 + k1/(4*pA) = 0.5 + 0.5/(4*0.1) = 0.5 + 1.25 = 1.75
  expected_lr <- 0.5 + 0.5 / (4 * 0.10)
  
  expect_equal(result$LR, expected_lr)
})

test_that("kinship_calculation correctly looks up full_siblings coefficients", {
  row <- create_genotype_row("D3S1358", "15", "16", "15", "16")  # AB-AB
  
  result <- kinship_calculation(
    row,
    mock_allele_freq,
    mock_kinship_matrix,
    tested_relationship = "full_siblings",
    tested_populations = "PopA"
  )
  
  # For AB-AB pattern with 2 shared alleles:
  # pA = 0.10, pB = 0.20
  # Rxp = 4*pA*pB/(pA+pB) = 4*0.1*0.2/0.3 = 0.08/0.3 = 0.2667
  # Rxu = 2*pA*pB = 0.04
  # LR = k0 + k1/Rxp + k2/Rxu = 0.25 + 0.5/0.2667 + 0.25/0.04
  pA <- 0.10
  pB <- 0.20
  Rxp <- (4 * pA * pB) / (pA + pB)
  Rxu <- 2 * pA * pB
  expected_lr <- 0.25 + 0.5 / Rxp + 0.25 / Rxu
  
  expect_equal(result$LR, expected_lr, tolerance = 1e-10)
  expect_equal(result$shared_alleles, 2)
  expect_equal(result$genotype_match, "AB-AB")
})

test_that("kinship_calculation correctly looks up cousins coefficients", {
  row <- create_genotype_row("D3S1358", "15", "15", "15", "16")  # AA-AB
  
  result <- kinship_calculation(
    row,
    mock_allele_freq,
    mock_kinship_matrix,
    tested_relationship = "cousins",
    tested_populations = "PopA"
  )
  
  # For AA-AB pattern with 1 shared allele:
  # LR = k0 + k1/(2*pA) = 0.875 + 0.125/(2*0.1) = 0.875 + 0.625 = 1.5
  expected_lr <- 0.875 + 0.125 / (2 * 0.10)
  
  expect_equal(result$LR, expected_lr)
  expect_equal(result$genotype_match, "AA-AB")
})

test_that("kinship_calculation correctly looks up unrelated coefficients", {
  row <- create_genotype_row("D3S1358", "15", "16", "15", "17")
  
  result <- kinship_calculation(
    row,
    mock_allele_freq,
    mock_kinship_matrix,
    tested_relationship = "unrelated",
    tested_populations = "PopA"
  )
  
  # For unrelated: k0=1, k1=0, k2=0
  # LR = 1 + 0 = 1 (regardless of genotype pattern)
  expect_equal(result$LR, 1)
})

# =============================================================================
# TEST: Allele frequency lookup
# =============================================================================

test_that("kinship_calculation uses correct population frequencies", {
  row <- create_genotype_row("D3S1358", "15", "16", "15", "17")
  
  # Test with PopA
  result_A <- kinship_calculation(
    row,
    mock_allele_freq,
    mock_kinship_matrix,
    tested_relationship = "half_siblings",
    tested_populations = "PopA"
  )
  
  # Test with PopB
  result_B <- kinship_calculation(
    row,
    mock_allele_freq,
    mock_kinship_matrix,
    tested_relationship = "half_siblings",
    tested_populations = "PopB"
  )
  
  # PopA: pA(15) = 0.10 → LR = 0.5 + 0.5/(4*0.1) = 1.75
  # PopB: pA(15) = 0.25 → LR = 0.5 + 0.5/(4*0.25) = 1.0
  expected_A <- 0.5 + 0.5 / (4 * 0.10)  # = 1.75
  expected_B <- 0.5 + 0.5 / (4 * 0.25)  # = 1.0
  
  expect_equal(result_A$LR, expected_A)
  expect_equal(result_B$LR, expected_B)
  expect_true(result_A$LR > result_B$LR)  # Rarer allele gives higher LR
})

test_that("kinship_calculation returns different LRs for different populations", {
  # Same genotype, same relationship, different frequency populations
  row <- create_genotype_row("D3S1358", "15", "16", "15", "16")  # AB-AB
  
  result_A <- kinship_calculation(
    row, mock_allele_freq, mock_kinship_matrix,
    "parent_child", "PopA"
  )
  
  result_B <- kinship_calculation(
    row, mock_allele_freq, mock_kinship_matrix,
    "parent_child", "PopB"
  )
  
  # LRs should be different because frequencies differ
  expect_false(result_A$LR == result_B$LR)
})

# =============================================================================
# TEST: Multiple tested populations
# =============================================================================

test_that("kinship_calculation handles multiple tested populations", {
  row <- create_genotype_row("D3S1358", "15", "16", "15", "17")
  
  result <- kinship_calculation(
    row,
    mock_allele_freq,
    mock_kinship_matrix,
    tested_relationship = "half_siblings",
    tested_populations = c("PopA", "PopB")
  )
  
  # Should return 2 rows, one for each population
  expect_equal(nrow(result), 2)
  expect_equal(result$tested_population, c("PopA", "PopB"))
  
  # Verify LRs are different
  expect_false(result$LR[1] == result$LR[2])
})

# =============================================================================
# TEST: All genotype patterns for half_siblings
# =============================================================================

test_that("half_siblings: 0 shared (AB-CD) gives LR = 0.5", {
  row <- create_genotype_row("D3S1358", "15", "16", "17", "14")  # Different alleles
  
  # Add allele 14 to mock frequencies
  mock_freq_extended <- rbind(
    mock_allele_freq,
    data.frame(population = "PopA", marker = "D3S1358", allele = "14", frequency = 0.40)
  )
  
  result <- kinship_calculation(
    row, mock_freq_extended, mock_kinship_matrix,
    "half_siblings", "PopA"
  )
  
  expect_equal(result$shared_alleles, 0)
  expect_equal(result$LR, 0.5)  # k0 for half_siblings
})

test_that("half_siblings: 1 shared AA-AB gives correct LR", {
  row <- create_genotype_row("D3S1358", "15", "15", "15", "16")  # AA vs AB
  
  result <- kinship_calculation(
    row, mock_allele_freq, mock_kinship_matrix,
    "half_siblings", "PopA"
  )
  
  # LR = k0 + k1/(2*pA) = 0.5 + 0.5/(2*0.1) = 0.5 + 2.5 = 3.0
  expected_lr <- 0.5 + 0.5 / (2 * 0.10)
  
  expect_equal(result$shared_alleles, 1)
  expect_equal(result$genotype_match, "AA-AB")
  expect_equal(result$LR, expected_lr)
})

test_that("half_siblings: 1 shared AB-AA gives correct LR", {
  row <- create_genotype_row("D3S1358", "15", "16", "15", "15")  # AB vs AA
  
  result <- kinship_calculation(
    row, mock_allele_freq, mock_kinship_matrix,
    "half_siblings", "PopA"
  )
  
  # LR = k0 + k1/(2*pA) = 0.5 + 0.5/(2*0.1) = 3.0
  expected_lr <- 0.5 + 0.5 / (2 * 0.10)
  
  expect_equal(result$shared_alleles, 1)
  expect_equal(result$genotype_match, "AB-AA")
  expect_equal(result$LR, expected_lr)
})

test_that("half_siblings: 1 shared AB-AC gives correct LR", {
  row <- create_genotype_row("D3S1358", "15", "16", "15", "17")  # AB vs AC
  
  result <- kinship_calculation(
    row, mock_allele_freq, mock_kinship_matrix,
    "half_siblings", "PopA"
  )
  
  # LR = k0 + k1/(4*pA) = 0.5 + 0.5/(4*0.1) = 0.5 + 1.25 = 1.75
  expected_lr <- 0.5 + 0.5 / (4 * 0.10)
  
  expect_equal(result$shared_alleles, 1)
  expect_equal(result$genotype_match, "AB-AC")
  expect_equal(result$LR, expected_lr)
})

test_that("half_siblings: 2 shared AA-AA gives correct LR (no k2 term)", {
  row <- create_genotype_row("D3S1358", "15", "15", "15", "15")  # AA vs AA
  
  result <- kinship_calculation(
    row, mock_allele_freq, mock_kinship_matrix,
    "half_siblings", "PopA"
  )
  
  # LR = k0 + k1/pA + k2/pA^2 = 0.5 + 0.5/0.1 + 0 = 0.5 + 5 = 5.5
  # Note: k2=0 for half_siblings, so no k2 term!
  expected_lr <- 0.5 + 0.5 / 0.10 + 0
  
  expect_equal(result$shared_alleles, 2)
  expect_equal(result$genotype_match, "AA-AA")
  expect_equal(result$LR, expected_lr)
})

test_that("half_siblings: 2 shared AB-AB gives correct LR (no k2 term)", {
  row <- create_genotype_row("D3S1358", "15", "16", "15", "16")  # AB vs AB
  
  result <- kinship_calculation(
    row, mock_allele_freq, mock_kinship_matrix,
    "half_siblings", "PopA"
  )
  
  # pA = 0.10, pB = 0.20
  # Rxp = 4*pA*pB/(pA+pB) = 0.08/0.3 = 0.2667
  # LR = k0 + k1/Rxp + k2/Rxu = 0.5 + 0.5/0.2667 + 0 = 2.375
  pA <- 0.10
  pB <- 0.20
  Rxp <- (4 * pA * pB) / (pA + pB)
  expected_lr <- 0.5 + 0.5 / Rxp + 0
  
  expect_equal(result$shared_alleles, 2)
  expect_equal(result$genotype_match, "AB-AB")
  expect_equal(result$LR, expected_lr, tolerance = 1e-10)
})

# =============================================================================
# TEST: Full siblings k2 term (differs from half_siblings)
# =============================================================================

test_that("full_siblings has higher LR than half_siblings for AA-AA due to k2 term", {
  row <- create_genotype_row("D3S1358", "15", "15", "15", "15")  # AA vs AA
  
  result_fs <- kinship_calculation(
    row, mock_allele_freq, mock_kinship_matrix,
    "full_siblings", "PopA"
  )
  
  result_hs <- kinship_calculation(
    row, mock_allele_freq, mock_kinship_matrix,
    "half_siblings", "PopA"
  )
  
  # FS: 0.25 + 0.5/0.1 + 0.25/0.01 = 0.25 + 5 + 25 = 30.25
  # HS: 0.5 + 0.5/0.1 + 0 = 5.5
  expect_gt(result_fs$LR, result_hs$LR)
  
  # Verify exact values
  expect_equal(result_fs$LR, 0.25 + 0.5/0.1 + 0.25/0.01)
  expect_equal(result_hs$LR, 0.5 + 0.5/0.1)
})

# =============================================================================
# TEST: Edge cases
# =============================================================================

test_that("kinship_calculation handles NA alleles", {
  row <- create_genotype_row("D3S1358", NA, "16", "15", "17")
  
  result <- kinship_calculation(
    row, mock_allele_freq, mock_kinship_matrix,
    "half_siblings", "PopA"
  )
  
  expect_true(is.na(result$LR))
  expect_true(is.na(result$shared_alleles))
})

test_that("kinship_calculation uses fallback frequency for unknown alleles", {
  # Allele "99" is not in the frequency table
  row <- create_genotype_row("D3S1358", "99", "99", "99", "99")  # AA vs AA with unknown allele
  
  result <- kinship_calculation(
    row, mock_allele_freq, mock_kinship_matrix,
    "half_siblings", "PopA"
  )
  
  # Should use FALLBACK_FREQ = 0.001
  # LR = 0.5 + 0.5/0.001 = 0.5 + 500 = 500.5
  expected_lr <- 0.5 + 0.5 / FALLBACK_FREQ
  
  expect_equal(result$LR, expected_lr)
})

test_that("kinship_calculation returns NA for missing locus in frequency data", {
  row <- create_genotype_row("UNKNOWN_LOCUS", "15", "16", "15", "17")
  
  result <- kinship_calculation(
    row, mock_allele_freq, mock_kinship_matrix,
    "half_siblings", "PopA"
  )
  
  expect_true(is.na(result$LR))
})

test_that("kinship_calculation errors on invalid relationship type", {
  row <- create_genotype_row("D3S1358", "15", "16", "15", "17")
  
  expect_error(
    kinship_calculation(
      row, mock_allele_freq, mock_kinship_matrix,
      "invalid_relationship", "PopA"
    ),
    "tested_relationship not found"
  )
})

# =============================================================================
# TEST: Cross-population comparison (key research question)
# =============================================================================

test_that("same pair gets different LRs with different population frequencies", {
  # This tests the core research question: using wrong population frequencies
  # affects the LR calculation
  
  row <- create_genotype_row("D3S1358", "15", "16", "15", "17",
                              population = "PopA", known_rel = "half_siblings")
  
  # Test with correct population (PopA)
  result_correct <- kinship_calculation(
    row, mock_allele_freq, mock_kinship_matrix,
    "half_siblings", "PopA"
  )
  
  # Test with wrong population (PopB)
  result_wrong <- kinship_calculation(
    row, mock_allele_freq, mock_kinship_matrix,
    "half_siblings", "PopB"
  )
  
  # LRs should differ because frequencies differ
  # PopA pA(15) = 0.10, PopB pA(15) = 0.25
  # Rarer allele in PopA gives higher LR
  expect_true(result_correct$LR != result_wrong$LR)
  expect_gt(result_correct$LR, result_wrong$LR)  # Rarer allele = higher LR
})

# =============================================================================
# TEST: Verify kinship coefficients match expected values
# =============================================================================

test_that("mock kinship matrix matches theoretical values", {
  # Parent-child: always share exactly 1 allele IBD
  pc <- mock_kinship_matrix |> filter(relationship_type == "parent_child")
  expect_equal(pc$k0, 0)
  expect_equal(pc$k1, 1)
  expect_equal(pc$k2, 0)
  
  # Full siblings: 25% share 0, 50% share 1, 25% share 2
  fs <- mock_kinship_matrix |> filter(relationship_type == "full_siblings")
  expect_equal(fs$k0, 0.25)
  expect_equal(fs$k1, 0.5)
  expect_equal(fs$k2, 0.25)
  
  # Half siblings (also avuncular, grandparent): 50% share 0, 50% share 1
  hs <- mock_kinship_matrix |> filter(relationship_type == "half_siblings")
  expect_equal(hs$k0, 0.5)
  expect_equal(hs$k1, 0.5)
  expect_equal(hs$k2, 0)
  
  # Cousins: 87.5% share 0, 12.5% share 1
  co <- mock_kinship_matrix |> filter(relationship_type == "cousins")
  expect_equal(co$k0, 0.875)
  expect_equal(co$k1, 0.125)
  expect_equal(co$k2, 0)
  
  # Unrelated: always share 0 IBD
  un <- mock_kinship_matrix |> filter(relationship_type == "unrelated")
  expect_equal(un$k0, 1)
  expect_equal(un$k1, 0)
  expect_equal(un$k2, 0)
})
