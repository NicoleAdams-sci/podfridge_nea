################################################################################
# Kinship‑LR utilities (long format version)
# Tina Lasisi, 2025‑06‑27
#
#  Required packages: dplyr, tibble, stringr, purrr (all tidyverse core)
#
#  Inputs (per function call)
#  --------------------------
#    row                    One tibble/data.frame row that contains
#                             ind1_allele1, ind1_allele2,
#                             ind2_allele1, ind2_allele2,
#                             locus
#
#    allele_frequency_data   Data frame with columns:
#                               population, marker, allele, frequency
#
#    kinship_matrix          Data frame with columns:
#                               relationship_type, k0, k1, k2
#
#    tested_relationship     Character scalar, e.g. "FullSibs"
#
#    tested_populations      Character vector, e.g. c("EAS", "EUR", "AFR")
#
#  Output
#  ------
#    A tibble with one row for each tested population, including
#      • all original row fields,
#      • shared_alleles                  (0, 1, or 2)
#      • genotype_match                  ("AA-AA", "AA-AB", "AB-AA", "AB-AB",
#                                         "AB-AC")
#      • tested_relationship             (copied input)
#      • tested_population               (one per row)
#      • LR                              (numeric, NA if indeterminable)
#
################################################################################

library(dplyr)
library(stringr)
library(purrr)

# -----------------------------------------------------------------------------#
# Global contingency allele frequency
# -----------------------------------------------------------------------------#
FALLBACK_FREQ <- 5 / (2 * 1036)   # ≈ 0.002414

# -----------------------------------------------------------------------------#
# Helper 1: Count shared alleles WITH multiplicity
# -----------------------------------------------------------------------------#
count_shared_alleles <- function(alleles1, alleles2) {
  # More robust type conversion
  alleles1 <- as.vector(as.character(alleles1))
  alleles2 <- as.vector(as.character(alleles2))
  
  # Remove any NA values
  alleles1 <- alleles1[!is.na(alleles1)]
  alleles2 <- alleles2[!is.na(alleles2)]
  
  if (length(alleles1) == 0 || length(alleles2) == 0) {
    return(0)
  }
  
  tbl1 <- table(alleles1)
  tbl2 <- table(alleles2)
  common <- intersect(names(tbl1), names(tbl2))
  
  if (length(common) == 0) {
    return(0)  # No shared alleles - this is the key fix!
  }
  
  sum(mapply(function(a) min(tbl1[[a]], tbl2[[a]]), common))
}

# -----------------------------------------------------------------------------#
# Helper 2: Build allele‑letter map and genotype strings
# -----------------------------------------------------------------------------#
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

# -----------------------------------------------------------------------------#
# Helper 3: LR formulae
# -----------------------------------------------------------------------------#
calculate_likelihood_ratio <- function(shared_alleles,
                                       genotype_match,
                                       pA = NA_real_, pB = NA_real_,
                                       k0, k1, k2) {
  
  # Probability that the tested relationship shares 0 IBD alleles
  if (shared_alleles == 0) return(k0)
  
  if (shared_alleles == 1) {
    Rxp <- switch(genotype_match,
                  "AA-AB" = 2 * pA,
                  "AB-AA" = 2 * pA,       # valid when ind1 is AB, ind2 is AA
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

# -----------------------------------------------------------------------------#
# Main driver (long format)
# -----------------------------------------------------------------------------#
kinship_calculation <- function(row,
                                allele_frequency_data,
                                kinship_matrix,
                                tested_relationship,
                                tested_populations) {
  
  # ----------------------------------------------------------------------
  # 0. Pre‑flight NA check
  # ----------------------------------------------------------------------
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
  
  # ----------------------------------------------------------------------
  # 1. Prepare allele character vectors
  # ----------------------------------------------------------------------
  alleles_ind1 <- as.character(c(row$ind1_allele1, row$ind1_allele2))
  alleles_ind2 <- as.character(c(row$ind2_allele1, row$ind2_allele2))
  
  # ----------------------------------------------------------------------
  # 2. Count shared alleles with multiplicity
  # ----------------------------------------------------------------------
  shared_alleles <- count_shared_alleles(alleles_ind1, alleles_ind2)
  
  # ----------------------------------------------------------------------
  # 3. Build deterministic labels and genotype string
  # ----------------------------------------------------------------------
  lab <- label_and_genotype(alleles_ind1, alleles_ind2)
  
  # ----------------------------------------------------------------------
  # 4. Retrieve k coefficients for the TESTED relationship
  # ----------------------------------------------------------------------
  kin_vals <- kinship_matrix |>
    filter(relationship_type == tested_relationship)
  
  if (nrow(kin_vals) != 1L)
    stop("tested_relationship not found or duplicated in kinship_matrix")
  
  k0 <- kin_vals$k0
  k1 <- kin_vals$k1
  k2 <- kin_vals$k2
  
  # ----------------------------------------------------------------------
  # 5. Helper to fetch allele frequencies with contingency
  # ----------------------------------------------------------------------
  fetch_freq <- function(freq_df, target_allele) {
    if (is.na(target_allele)) return(NA_real_)
    idx <- match(target_allele, freq_df$allele)
    if (is.na(idx)) return(FALLBACK_FREQ)         # allele not present
    f <- freq_df$frequency[[idx]]
    if (is.na(f) || f == 0) FALLBACK_FREQ else f
  }
  
  # ----------------------------------------------------------------------
  # 6. Iterate over tested populations → build long‑form rows
  # ----------------------------------------------------------------------
  rows_out <- map_dfr(tested_populations, function(pop) {
    
    pop_freqs <- allele_frequency_data |>
      filter(population == pop, marker == row$locus)
    
    # If locus missing entirely, LR NA
    if (nrow(pop_freqs) == 0L) {
      LR_val <- NA_real_
    } else {
      # Map letters → allele strings
      A_val <- if ("A" %in% names(lab$allele_map)) lab$allele_map[["A"]] else NA_character_
      B_val <- if ("B" %in% names(lab$allele_map)) lab$allele_map[["B"]] else NA_character_
      
      # Fetch freqs (with contingency)
      pA <- fetch_freq(pop_freqs, A_val)
      pB <- fetch_freq(pop_freqs, B_val)
      
      # If pA could not be determined, LR is NA
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
    
    # Assemble one output row for this population
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

################################################################################
# End of file
################################################################################