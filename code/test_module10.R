##### Module 10 Standalone Test #####
# USAGE: Rscript code/test_module10.R
#
# Purpose:
#   Tests module10_ranking_db_assembler.R in isolation by:
#     1. Generating a tiny focal family (module 7) — 2 focal individuals
#     2. Generating a tiny unrelated pool (module 8) — 10 individuals per pop
#     3. Running assemble_ranking_database() for a single focal + relative
#     4. Running assemble_ranking_replicates() across all focal individuals
#     5. Checking output structure, dimensions, and key column values
#
# Expected output:
#   - Assembled paired data frame with (10*4 + 1) = 41 database members
#   - One row per database member per locus (41 * 29 = 1189 rows)
#   - Columns: batch_id, pair_id, individual_id, population, locus,
#              focal_A1, focal_A2, ind2_A1, ind2_A2,
#              known_relationship, is_true_relative
#   - Exactly 1 row per locus where is_true_relative == TRUE
#   - population column is "Asian" for all rows

cat("=== Module 10 Standalone Test ===\n")
cat("Started at:", format(Sys.time()), "\n\n")

# ------------------------------------------------------------------------------
# 0. Load packages and source modules
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(purrr)
})

source("code/LR_kinship_utility_functions.R")
source("code/module1_allele_simulator.R")
source("code/module2_STR_profile_simulator.R")
source("code/module3_related_individual_simulator.R")
source("code/module7_single_pop_focal_generator.R")
source("code/module8_unrelated_pool_generator.R")
source("code/module10_ranking_db_assembler.R")

# Note: must be named kinship_matrix — module 3 references this as a global variable
kinship_matrix <- fread("data/kinship_coefficients.csv")

# Small test parameters
N_FOCAL          <- 2
N_UNRELATED      <- 10   # per population, so 40 total
FOCAL_POPULATION <- "Asian"
RELATIVE_TYPE    <- "parent_child"

cat(sprintf("Test parameters:\n"))
cat(sprintf("  Focal individuals    : %d\n", N_FOCAL))
cat(sprintf("  Unrelated per pop    : %d (total: %d)\n", N_UNRELATED, N_UNRELATED * 4))
cat(sprintf("  Focal population     : %s\n", FOCAL_POPULATION))
cat(sprintf("  Relative type tested : %s\n\n", RELATIVE_TYPE))

# ------------------------------------------------------------------------------
# STEP 1: Generate tiny focal family (module 7)
# ------------------------------------------------------------------------------

cat("--- STEP 1: Generate focal family (module 7) ---\n")

focal_result <- generate_single_pop_focal(
  population            = FOCAL_POPULATION,
  n_focal               = N_FOCAL,
  relationship_counts   = list(parent_child = 1, full_siblings = 1),
  loci_list             = loci_list,
  allele_frequency_data = df_allelefreq,
  kinship_coefficients  = kinship_matrix,
  output_dir            = "output/focal_ranking_test"
)

focal_family_data <- focal_result$data

cat(sprintf("Rows in focal_family_data     : %d\n", nrow(focal_family_data)))
cat(sprintf("Unique focal_ids              : %s\n",
            paste(unique(focal_family_data$focal_id), collapse = ", ")))
cat(sprintf("Unique relationship_to_focal  : %s\n",
            paste(unique(focal_family_data$relationship_to_focal), collapse = ", ")))
cat(sprintf("Unique loci                   : %d\n\n",
            length(unique(focal_family_data$locus))))

# ------------------------------------------------------------------------------
# STEP 2: Generate tiny unrelated pool (module 8)
# ------------------------------------------------------------------------------

cat("--- STEP 2: Generate unrelated pool (module 8) ---\n")

unrelated_results <- generate_multiple_pop_unrelated(
  populations           = c("AfAm", "Cauc", "Hispanic", "Asian"),
  n_unrelated_per_pop   = N_UNRELATED,
  loci_list             = loci_list,
  allele_frequency_data = df_allelefreq,
  output_dir            = "output/unrelated_pool",
  use_single_datetime   = TRUE
)

unrelated_pool_data <- map_dfr(unrelated_results$file_path, fread)

# Note: individual_ids repeat across populations in module 8 output (e.g. unrel_001)
# so count by nrow / n_loci to get true individual count
n_loci_check <- length(unique(unrelated_pool_data$locus))
cat(sprintf("Total unrelated individuals   : %d\n",
            nrow(unrelated_pool_data) / n_loci_check))
cat(sprintf("Populations in pool           : %s\n\n",
            paste(unique(unrelated_pool_data$population), collapse = ", ")))

# ------------------------------------------------------------------------------
# STEP 3: Test assemble_ranking_database() — single focal, single relative
# ------------------------------------------------------------------------------

cat("--- STEP 3: Test assemble_ranking_database() (single replicate) ---\n")

first_focal_id <- unique(focal_family_data$focal_id[
  focal_family_data$relationship_to_focal == "self"
])[1]

cat(sprintf("Using focal_id: %s\n", first_focal_id))

paired_db <- assemble_ranking_database(
  focal_family_data   = focal_family_data,
  unrelated_pool_data = unrelated_pool_data,
  focal_id            = first_focal_id,
  relative_type       = RELATIVE_TYPE,
  relative_index      = 1
)

# --- Checks ---
n_loci        <- length(unique(paired_db$locus))
n_db_members  <- length(unique(paired_db$individual_id))
n_true_rel    <- length(unique(paired_db$individual_id[paired_db$is_true_relative == TRUE]))
expected_rows <- n_db_members * n_loci

cat(sprintf("\nOutput checks:\n"))
cat(sprintf("  Columns present       : %s\n", paste(names(paired_db), collapse = ", ")))
cat(sprintf("  Total rows            : %d (expected: %d)\n", nrow(paired_db), expected_rows))
cat(sprintf("  Unique database members : %d (expected: %d)\n",
            n_db_members, N_UNRELATED * 4 + 1))
cat(sprintf("  Unique loci           : %d\n", n_loci))
cat(sprintf("  True relatives        : %d (expected: 1)\n", n_true_rel))
cat(sprintf("  population values     : %s (expected: all %s)\n",
            paste(unique(paired_db$population), collapse = ", "),
            FOCAL_POPULATION))
cat(sprintf("  known_relationship values : %s\n",
            paste(unique(paired_db$known_relationship), collapse = ", ")))

# Check no NA alleles
n_na_focal <- sum(is.na(paired_db$focal_A1) | is.na(paired_db$focal_A2))
n_na_ind2  <- sum(is.na(paired_db$ind2_A1)  | is.na(paired_db$ind2_A2))
cat(sprintf("  NA focal alleles      : %d (expected: 0)\n", n_na_focal))
cat(sprintf("  NA ind2 alleles       : %d (expected: 0)\n", n_na_ind2))

# Peek at first few rows
cat("\nFirst 5 rows of paired_db:\n")
print(head(paired_db, 5))

# Peek at the true relative rows
cat("\nTrue relative rows (first locus only):\n")
print(paired_db[paired_db$is_true_relative == TRUE, ][1, ])

# ------------------------------------------------------------------------------
# STEP 4: Test assemble_ranking_replicates() — all focal individuals
# ------------------------------------------------------------------------------

cat("\n--- STEP 4: Test assemble_ranking_replicates() (all replicates) ---\n")

replicate_dbs <- assemble_ranking_replicates(
  focal_family_data   = focal_family_data,
  unrelated_pool_data = unrelated_pool_data,
  relative_type       = RELATIVE_TYPE,
  relative_index      = 1
)

cat(sprintf("\nOutput checks:\n"))
cat(sprintf("  Number of replicates  : %d (expected: %d)\n",
            length(replicate_dbs), N_FOCAL))
cat(sprintf("  Replicate names       : %s\n",
            paste(names(replicate_dbs), collapse = ", ")))
cat(sprintf("  Rows per replicate    : %s\n",
            paste(sapply(replicate_dbs, nrow), collapse = ", ")))

# Check each replicate has exactly one true relative
n_true_per_rep <- sapply(replicate_dbs, function(db) {
  length(unique(db$individual_id[db$is_true_relative == TRUE]))
})
cat(sprintf("  True relatives per replicate : %s (expected: all 1)\n",
            paste(n_true_per_rep, collapse = ", ")))

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------

cat("\n=== Module 10 Test Summary ===\n")

checks <- c(
  "Correct number of rows"      = nrow(paired_db) == expected_rows,
  "Correct database size"       = n_db_members == (N_UNRELATED * 4 + 1),
  "Exactly one true relative"   = n_true_rel == 1,
  "Population label consistent" = all(paired_db$population == FOCAL_POPULATION),
  "No NA focal alleles"         = n_na_focal == 0,
  "No NA ind2 alleles"          = n_na_ind2 == 0,
  "Correct replicate count"     = length(replicate_dbs) == N_FOCAL,
  "One true relative per rep"   = all(n_true_per_rep == 1)
)

for (check_name in names(checks)) {
  status <- ifelse(checks[check_name], "PASS", "FAIL")
  cat(sprintf("  [%s] %s\n", status, check_name))
}

if (all(checks)) {
  cat("\nSUCCESS: All checks passed. Module 10 is working correctly.\n")
} else {
  cat("\nWARNING: Some checks failed. Review output above.\n")
}

cat("Completed at:", format(Sys.time()), "\n")
