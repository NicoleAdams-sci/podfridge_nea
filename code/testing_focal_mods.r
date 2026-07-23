##### Compare calculate_single_locus_lr() vs calculate_single_locus_lr_fast() #####
# Run this interactively or via Rscript from the project root
# (/nfs/turbo/lsa-tlasisi1/podfridge_simulations/)

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(tidyverse)
})

setwd("/nfs/turbo/lsa-tlasisi1/podfridge_simulations/")

source("code/module11_ranking_lr_calculator_fast.R")   # brings in df_allelefreq via module4 source chain
source("code/module12_ranking_outcome_recorder.R")
source("code/focal_test_helper_fns.R")
source("code/module4_single_locus_LR_fast.R")      # the fast function, in case module11 wasn't swapped

kinship_matrix <- fread("data/kinship_coefficients.csv")

core_loci <- fread("data/core_CODIS_loci.csv")
loci_lists <- lapply(c("core_13", "expanded_20"), function(col) {
  core_loci |> filter(get(col) == 1) |> pull(locus)
})
names(loci_lists) <- c("core_13", "expanded_20")
loci_needed <- unique(unlist(loci_lists))

# -------------------------------------------------------------------------
# Same files as your c10_001 test run
# -------------------------------------------------------------------------
PAIR_FILE           <- "output/pairs_all_parent_child_n1000_chunk10_20260130.csv"
UNRELATED_POOL_FILE <- "output/unrelated_pool/all_N10000_combined_unrelated_20260702_153012.csv"

TARGET_BATCH_ID <- "20260130_135230"
TARGET_PAIR_ID  <- "c10_001"

pair_data <- fread(
  PAIR_FILE,
  select = c("batch_id", "pair_id", "population", "known_relationship",
             "locus", "focal_A1", "focal_A2")
) |> filter(locus %in% loci_needed)

unrelated_pool_data <- fread(
  UNRELATED_POOL_FILE,
  select = c("batch_id", "database_composition", "database_label", "individual_id",
             "relationship_to_focal", "source_frequency_population",
             "population", "locus", "A1", "A2")
) |> filter(locus %in% loci_needed)

pair_data_one_pair <- pair_data |>
  filter(batch_id == TARGET_BATCH_ID, pair_id == TARGET_PAIR_ID)

# -------------------------------------------------------------------------
# Build the SAME paired_unrelated_db the real run computes, then slice small
# -------------------------------------------------------------------------
paired_unrelated_db <- assemble_unrelated_database_from_existing_pair(
  pair_data_one_pair  = pair_data_one_pair,
  unrelated_pool_data = unrelated_pool_data,
  batch_id            = paste0("pair_", TARGET_PAIR_ID, "_unrelateds")
)

small_paired_db <- paired_unrelated_db |>
  filter(pair_id %in% head(unique(pair_id), 20))

cat(sprintf("small_paired_db: %d rows, %d candidates, %d loci\n",
            nrow(small_paired_db), n_distinct(small_paired_db$pair_id),
            n_distinct(small_paired_db$locus)))

# -------------------------------------------------------------------------
# Compare old vs fast for full_siblings (the relationship showing divergence)
# -------------------------------------------------------------------------
old_lr <- calculate_single_locus_lr(
  pair_data = small_paired_db,
  tested_relationship = "full_siblings",
  tested_populations = "all",
  allele_frequency_data = df_allelefreq,
  kinship_coefficients = kinship_matrix
)

new_lr <- calculate_single_locus_lr_fast(
  pair_data = small_paired_db,
  tested_relationship = "full_siblings",
  tested_populations = "all",
  allele_frequency_data = df_allelefreq,
  kinship_coefficients = kinship_matrix
)

compare <- old_lr |>
  select(pair_id, locus, LR_old = LR) |>
  left_join(new_lr |> select(pair_id, locus, LR_new = LR), by = c("pair_id", "locus")) |>
  mutate(diff = abs(LR_old - LR_new))

cat("\n--- Rows where old and fast disagree (full_siblings) ---\n")
print(compare |> filter(diff > 1e-9 | is.na(diff)))

cat(sprintf("\nTotal rows: %d | Mismatches: %d\n",
            nrow(compare), sum(compare$diff > 1e-9 | is.na(compare$diff))))

# -------------------------------------------------------------------------
# Also check parent_child for completeness
# -------------------------------------------------------------------------
old_lr_pc <- calculate_single_locus_lr(
  pair_data = small_paired_db,
  tested_relationship = "parent_child",
  tested_populations = "all",
  allele_frequency_data = df_allelefreq,
  kinship_coefficients = kinship_matrix
)

new_lr_pc <- calculate_single_locus_lr_fast(
  pair_data = small_paired_db,
  tested_relationship = "parent_child",
  tested_populations = "all",
  allele_frequency_data = df_allelefreq,
  kinship_coefficients = kinship_matrix
)

compare_pc <- old_lr_pc |>
  select(pair_id, locus, LR_old = LR) |>
  left_join(new_lr_pc |> select(pair_id, locus, LR_new = LR), by = c("pair_id", "locus")) |>
  mutate(diff = abs(LR_old - LR_new))

cat("\n--- Rows where old and fast disagree (parent_child) ---\n")
print(compare_pc |> filter(diff > 1e-9 | is.na(diff)))

cat(sprintf("\nTotal rows: %d | Mismatches: %d\n",
            nrow(compare_pc), sum(compare_pc$diff > 1e-9 | is.na(compare_pc$diff))))