library(tidyverse)
library(data.table)

setwd("/nfs/turbo/lsa-tlasisi1/podfridge_simulations/")

source("code/module11_ranking_lr_calculator_fast.R")
source("code/module12_ranking_outcome_recorder.R")
source("code/focal_test_helper_fns.R")

kinship_matrix <- fread("data/kinship_coefficients.csv")

tested_relationships <- c("parent_child", "full_siblings")
tested_populations <- "all"

core_loci <- fread("data/core_CODIS_loci.csv")
columns <- c("core_13", "identifiler_15", "expanded_20", "supplementary")
loci_lists <- lapply(columns, function(col) {
  core_loci |>
    filter(get(col) == 1) |>
    pull(locus)
})
names(loci_lists) <- columns
loci_sets_to_use <- c("core_13", "expanded_20")
loci_sets <- list(
  core_13 = loci_lists$core_13,
  expanded_20 = loci_lists$expanded_20
)

# Extract all unique loci from allele frequency data
loci_list <- df_allelefreq |> pull(marker) |> unique()
loci_lists$autosomal_29 <- loci_list

#PAIR_FILE <- "output/pairs_all_parent_child_n1000_chunk3_20260130.csv"
#COMBINED_TRUE_LR_FILE <- "output/combined_LR/combined_LR_all_parent_child_n1000_chunk3_20260130.csv"
PAIR_FILE <- "output/pairs_all_cousins_n1000_chunk3_20260130.csv"
COMBINED_TRUE_LR_FILE <- "output/combined_LR/combined_LR_all_cousins_n1000_chunk3_20260130.csv"
UNRELATED_POOL_FILE <- "output/unrelated_pool/all_N1000_combined_unrelated_20260702_112423.csv"

pair_data <- fread(PAIR_FILE)
combined_true_lr_data <- fread(COMBINED_TRUE_LR_FILE)
unrelated_pool_data <- fread(UNRELATED_POOL_FILE)

loci_needed <- unique(unlist(loci_sets))

pair_data <- pair_data |>
  filter(locus %in% loci_needed)

unrelated_pool_data <- unrelated_pool_data |>
  filter(locus %in% loci_needed)

TARGET_PAIR_ID <- "c03_001"

#Subset to one pair
pair_data_one_pair <- pair_data |>
  filter(pair_id == TARGET_PAIR_ID)

combined_true_lr_one_pair <- combined_true_lr_data |>
  filter(pair_id == TARGET_PAIR_ID)

true_combined <- format_existing_true_combined_lr(
  combined_lr_one_pair = combined_true_lr_one_pair,
  tested_relationships = tested_relationships,
  tested_populations = tested_populations,
  loci_sets_to_use = loci_sets_to_use
)

paired_unrelated_db <- assemble_unrelated_database_from_existing_pair(
  pair_data_one_pair = pair_data_one_pair,
  unrelated_pool_data = unrelated_pool_data,
  batch_id = paste0("pair_", TARGET_PAIR_ID, "_unrelateds")
)

unrelated_combined <- calculate_ranking_lrs(
  paired_db = paired_unrelated_db,
  tested_relationships = tested_relationships,
  tested_populations = tested_populations,
  loci_sets = loci_sets,
  allele_frequency_data = df_allelefreq,
  kinship_coefficients = kinship_matrix,
  output_dir = "output/focal_ranking_test",
  save_results = FALSE
)

unrelated_combined$focal_id <- TARGET_PAIR_ID

ranking_input <- bind_rows(
  unrelated_combined,
  true_combined
)

ranking <- record_ranking_outcomes(
  lr_results = ranking_input,
  top_n = 200,
  output_dir = "output/focal_ranking_test",
  save_results = FALSE
)

outcomes <- ranking$outcomes |>
  mutate(
    original_pair_id = TARGET_PAIR_ID,
    true_relationship = known_relationship,
    top_200 = rank <= 200,
    top_100 = rank <= 100,
    top_50  = rank <= 50,
    top_10  = rank <= 10
  ) |>
  select(
    original_pair_id,
    true_relationship,
    loci_set,
    tested_relationship,
    tested_population,
    individual_id,
    combined_LR,
    rank,
    n_database,
    top_200,
    top_100,
    top_50,
    top_10
  )
