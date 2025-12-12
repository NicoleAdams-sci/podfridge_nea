# Test kinship_calculation()
# This function combines the three pure functions and looks up frequencies/kinship coefficients

library(testthat)
library(dplyr)
library(purrr)

# =============================================================================
# Define helper functions locally (copied from LR_kinship_utility_functions.R)
# =============================================================================

FALLBACK_FREQ <- 5 / (2 * 1036)  # ≈ 0.002414

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
  shared_vals      <- intersect(alleles1, alleles2)
  unique_ind1_vals <- setdiff(alleles1, shared_vals)
  unique_ind2_vals <- setdiff(alleles2, c(shared_vals, unique_ind1_vals))
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
  allele_map <- unlist(allele_map)
  translate <- function(vec) {
    letters <- names(allele_map)[match(vec, allele_map)]
    paste(sort(letters), collapse = "")
  }
  geno1 <- translate(alleles1)
  geno2 <- translate(alleles2)
  genotype_match <- paste(geno1, geno2, sep = "-")
  list(allele_map = allele_map, genotype_ind1 = geno1, genotype_ind2 = geno2, genotype_match = genotype_match)
}

calculate_likelihood_ratio <- function(shared_alleles, genotype_match, pA = NA_real_, pB = NA_real_, k0, k1, k2) {
  if (shared_alleles == 0) return(k0)
  if (shared_alleles == 1) {
    Rxp <- switch(genotype_match,
                  "AA-AB" = 2 * pA,
                  "AB-AA" = 2 * pA,
                  "AB-AC" = 4 * pA,
                  stop("Invalid genotype_match for 1 shared allele."))
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

# Main function to test
kinship_calculation <- function(row, allele_frequency_data, kinship_matrix, tested_relationship, tested_populations) {
  if (anyNA(row[c("ind1_allele1", "ind1_allele2", "ind2_allele1", "ind2_allele2")])) {
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
  
  kin_vals <- kinship_matrix |> filter(relationship_type == tested_relationship)
  if (nrow(kin_vals) != 1L) stop("tested_relationship not found or duplicated in kinship_matrix")
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
    pop_freqs <- allele_frequency_data |> filter(population == pop, marker == row$locus)
    if (nrow(pop_freqs) == 0L) {
      LR_val <- NA_real_
    } else {
      A_val <- if ("A" %in% names(lab$allele_map)) lab$allele_map[["A"]] else NA_character_
      B_val <- if ("B" %in% names(lab$allele_map)) lab$allele_map[["B"]] else NA_character_
      pA <- fetch_freq(pop_freqs, A_val)
      pB <- fetch_freq(pop_freqs, B_val)
      if (is.na(pA)) {
        LR_val <- NA_real_
      } else if (shared_alleles == 2 && lab$genotype_match == "AB-AB" && is.na(pB)) {
        LR_val <- NA_real_
      } else {
        LR_val <- calculate_likelihood_ratio(shared_alleles, lab$genotype_match, pA, pB, k0, k1, k2)
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
# MOCK DATA SETUP
# =============================================================================

# Mock allele frequency data - simple controlled values
mock_allele_freq <- data.frame(
  population = c("PopA", "PopA", "PopA", "PopA",
                 "PopB", "PopB", "PopB", "PopB"),
  marker = c("LOCUS1", "LOCUS1", "LOCUS1", "LOCUS1",
             "LOCUS1", "LOCUS1", "LOCUS1", "LOCUS1"),
  allele = c("15", "16", "17", "18",
             "15", "16", "17", "18"),
  frequency = c(0.25, 0.30, 0.20, 0.25,   # PopA frequencies
                0.10, 0.40, 0.30, 0.20),  # PopB frequencies (different!)
  stringsAsFactors = FALSE
)

# Mock kinship coefficients
mock_kinship <- data.frame(
  relationship_type = c("parent_child", "full_siblings", "unrelated"),
  k0 = c(0, 0.25, 1),
  k1 = c(1, 0.5, 0),
  k2 = c(0, 0.25, 0),
  stringsAsFactors = FALSE
)

# =============================================================================
# TEST: Basic LR calculation with known frequencies
# =============================================================================

test_that("kinship_calculation returns correct LR for AB-AB parent-child with PopA frequencies", {
  # Genotype: (15, 16) vs (15, 16) = AB-AB, 2 shared
  # PopA: pA(15)=0.25, pB(16)=0.30
  # Parent-child: k0=0, k1=1, k2=0
  # Formula: LR = k0 + k1*(pA+pB)/(4*pA*pB) + k2/(2*pA*pB)
  #            = 0 + 1*(0.55)/(0.30) + 0 = 1.8333...
  
  row <- list(
    ind1_allele1 = "15", ind1_allele2 = "16",
    ind2_allele1 = "15", ind2_allele2 = "16",
    locus = "LOCUS1"
  )
  
  result <- kinship_calculation(
    row = row,
    allele_frequency_data = mock_allele_freq,
    kinship_matrix = mock_kinship,
    tested_relationship = "parent_child",
    tested_populations = "PopA"
  )
  
  expected_lr <- (0.25 + 0.30) / (4 * 0.25 * 0.30)  # = 1.8333...
  
  expect_equal(nrow(result), 1)
  expect_equal(result$tested_population, "PopA")
  expect_equal(result$tested_relationship, "parent_child")
  expect_equal(result$shared_alleles, 2)
  expect_equal(result$genotype_match, "AB-AB")
  expect_equal(result$LR, expected_lr, tolerance = 1e-10)
})

test_that("same genotype gives different LR with different population frequencies", {
  # Same genotype (15, 16) vs (15, 16)
  # But PopA has pA=0.25, pB=0.30 and PopB has pA=0.10, pB=0.40
  # LR should be different!
  
  row <- list(
    ind1_allele1 = "15", ind1_allele2 = "16",
    ind2_allele1 = "15", ind2_allele2 = "16",
    locus = "LOCUS1"
  )
  
  result <- kinship_calculation(
    row = row,
    allele_frequency_data = mock_allele_freq,
    kinship_matrix = mock_kinship,
    tested_relationship = "parent_child",
    tested_populations = c("PopA", "PopB")
  )
  
  expect_equal(nrow(result), 2)
  
  lr_popA <- result$LR[result$tested_population == "PopA"]
  lr_popB <- result$LR[result$tested_population == "PopB"]
  
  # Hand-calculate expected values
  expected_lr_popA <- (0.25 + 0.30) / (4 * 0.25 * 0.30)  # pA=0.25, pB=0.30
  expected_lr_popB <- (0.10 + 0.40) / (4 * 0.10 * 0.40)  # pA=0.10, pB=0.40
  
  expect_equal(lr_popA, expected_lr_popA, tolerance = 1e-10)
  expect_equal(lr_popB, expected_lr_popB, tolerance = 1e-10)
  expect_false(lr_popA == lr_popB)  # They should be different
})

# =============================================================================
# TEST: AA-AA pattern
# =============================================================================

test_that("kinship_calculation returns correct LR for AA-AA parent-child", {
  # Genotype: (15, 15) vs (15, 15) = AA-AA, 2 shared
  # PopA: pA(15)=0.25
  # Parent-child: k0=0, k1=1, k2=0
  # Formula: LR = k0 + k1/pA + k2/pA^2 = 0 + 1/0.25 + 0 = 4.0
  
  row <- list(
    ind1_allele1 = "15", ind1_allele2 = "15",
    ind2_allele1 = "15", ind2_allele2 = "15",
    locus = "LOCUS1"
  )
  
  result <- kinship_calculation(
    row = row,
    allele_frequency_data = mock_allele_freq,
    kinship_matrix = mock_kinship,
    tested_relationship = "parent_child",
    tested_populations = "PopA"
  )
  
  expect_equal(result$genotype_match, "AA-AA")
  expect_equal(result$LR, 4.0)  # 1/0.25
})

test_that("kinship_calculation returns correct LR for AA-AA full siblings", {
  # Full siblings: k0=0.25, k1=0.5, k2=0.25
  # pA=0.25
  # LR = 0.25 + 0.5/0.25 + 0.25/0.0625 = 0.25 + 2 + 4 = 6.25
  
  row <- list(
    ind1_allele1 = "15", ind1_allele2 = "15",
    ind2_allele1 = "15", ind2_allele2 = "15",
    locus = "LOCUS1"
  )
  
  result <- kinship_calculation(
    row = row,
    allele_frequency_data = mock_allele_freq,
    kinship_matrix = mock_kinship,
    tested_relationship = "full_siblings",
    tested_populations = "PopA"
  )
  
  expect_equal(result$LR, 6.25)
})

# =============================================================================
# TEST: AA-AB pattern (1 shared allele)
# =============================================================================

test_that("kinship_calculation returns correct LR for AA-AB parent-child", {
  # Genotype: (15, 15) vs (15, 16) = AA-AB, 1 shared
  # PopA: pA(15)=0.25
  # Parent-child: k0=0, k1=1, k2=0
  # Formula: LR = k0 + k1/(2*pA) = 0 + 1/(2*0.25) = 2.0
  
  row <- list(
    ind1_allele1 = "15", ind1_allele2 = "15",
    ind2_allele1 = "15", ind2_allele2 = "16",
    locus = "LOCUS1"
  )
  
  result <- kinship_calculation(
    row = row,
    allele_frequency_data = mock_allele_freq,
    kinship_matrix = mock_kinship,
    tested_relationship = "parent_child",
    tested_populations = "PopA"
  )
  
  expect_equal(result$genotype_match, "AA-AB")
  expect_equal(result$shared_alleles, 1)
  expect_equal(result$LR, 2.0)
})

# =============================================================================
# TEST: AB-AC pattern (1 shared allele)
# =============================================================================

test_that("kinship_calculation returns correct LR for AB-AC parent-child", {
  # Genotype: (15, 16) vs (15, 17) = AB-AC, 1 shared
  # PopA: pA(15)=0.25
  # Parent-child: k0=0, k1=1, k2=0
  # Formula: LR = k0 + k1/(4*pA) = 0 + 1/(4*0.25) = 1.0
  
  row <- list(
    ind1_allele1 = "15", ind1_allele2 = "16",
    ind2_allele1 = "15", ind2_allele2 = "17",
    locus = "LOCUS1"
  )
  
  result <- kinship_calculation(
    row = row,
    allele_frequency_data = mock_allele_freq,
    kinship_matrix = mock_kinship,
    tested_relationship = "parent_child",
    tested_populations = "PopA"
  )
  
  expect_equal(result$genotype_match, "AB-AC")
  expect_equal(result$shared_alleles, 1)
  expect_equal(result$LR, 1.0)
})

# =============================================================================
# TEST: 0 shared alleles
# =============================================================================

test_that("kinship_calculation returns k0 for 0 shared alleles", {
  # Genotype: (15, 16) vs (17, 18) = AB-CD, 0 shared
  # Parent-child: k0=0 -> LR = 0 (exclusion)
  # Unrelated: k0=1 -> LR = 1
  
  row <- list(
    ind1_allele1 = "15", ind1_allele2 = "16",
    ind2_allele1 = "17", ind2_allele2 = "18",
    locus = "LOCUS1"
  )
  
  result_pc <- kinship_calculation(
    row = row,
    allele_frequency_data = mock_allele_freq,
    kinship_matrix = mock_kinship,
    tested_relationship = "parent_child",
    tested_populations = "PopA"
  )
  
  result_unrel <- kinship_calculation(
    row = row,
    allele_frequency_data = mock_allele_freq,
    kinship_matrix = mock_kinship,
    tested_relationship = "unrelated",
    tested_populations = "PopA"
  )
  
  expect_equal(result_pc$shared_alleles, 0)
  expect_equal(result_pc$LR, 0)  # Parent-child exclusion
  expect_equal(result_unrel$LR, 1)  # Unrelated
})

# =============================================================================
# TEST: Fallback frequency for unknown alleles
# =============================================================================

test_that("kinship_calculation uses FALLBACK_FREQ for unknown alleles", {
  # Allele "99" is not in our mock frequency data
  # Should use FALLBACK_FREQ ≈ 0.002414
  
  row <- list(
    ind1_allele1 = "99", ind1_allele2 = "99",
    ind2_allele1 = "99", ind2_allele2 = "99",
    locus = "LOCUS1"
  )
  
  result <- kinship_calculation(
    row = row,
    allele_frequency_data = mock_allele_freq,
    kinship_matrix = mock_kinship,
    tested_relationship = "parent_child",
    tested_populations = "PopA"
  )
  
  # AA-AA with pA = FALLBACK_FREQ
  # LR = 1 / FALLBACK_FREQ
  expected_lr <- 1 / FALLBACK_FREQ
  
  expect_equal(result$genotype_match, "AA-AA")
  expect_equal(result$LR, expected_lr, tolerance = 1e-10)
})

# =============================================================================
# TEST: NA handling
# =============================================================================

test_that("kinship_calculation returns NA for rows with NA alleles", {
  row <- list(
    ind1_allele1 = "15", ind1_allele2 = NA,
    ind2_allele1 = "15", ind2_allele2 = "16",
    locus = "LOCUS1"
  )
  
  result <- kinship_calculation(
    row = row,
    allele_frequency_data = mock_allele_freq,
    kinship_matrix = mock_kinship,
    tested_relationship = "parent_child",
    tested_populations = "PopA"
  )
  
  expect_true(is.na(result$LR))
  expect_true(is.na(result$shared_alleles))
})

# =============================================================================
# TEST: Multiple populations at once
# =============================================================================

test_that("kinship_calculation returns one row per tested population", {
  row <- list(
    ind1_allele1 = "15", ind1_allele2 = "16",
    ind2_allele1 = "15", ind2_allele2 = "16",
    locus = "LOCUS1"
  )
  
  result <- kinship_calculation(
    row = row,
    allele_frequency_data = mock_allele_freq,
    kinship_matrix = mock_kinship,
    tested_relationship = "parent_child",
    tested_populations = c("PopA", "PopB")
  )
  
  expect_equal(nrow(result), 2)
  expect_setequal(result$tested_population, c("PopA", "PopB"))
})

# =============================================================================
# TEST: Missing locus returns NA
# =============================================================================

test_that("kinship_calculation returns NA LR for missing locus", {
  row <- list(
    ind1_allele1 = "15", ind1_allele2 = "16",
    ind2_allele1 = "15", ind2_allele2 = "16",
    locus = "NONEXISTENT_LOCUS"
  )
  
  result <- kinship_calculation(
    row = row,
    allele_frequency_data = mock_allele_freq,
    kinship_matrix = mock_kinship,
    tested_relationship = "parent_child",
    tested_populations = "PopA"
  )
  
  expect_true(is.na(result$LR))
})
