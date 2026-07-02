##### Run Ranking for One Chunk of Existing Pair IDs #####
#
# Usage:
#   Rscript code/run_pair_ranking_chunk.R \
#     <PAIR_FILE> \
#     <COMBINED_TRUE_LR_FILE> \
#     <UNRELATED_POOL_FILE> \
#     <CHUNK_INDEX> \
#     <CHUNK_SIZE> \
#     <OUTPUT_DIR>
#
# Example:
#   Rscript code/run_pair_ranking_chunk.R \
#     output/pairs_all_cousins_n1000_chunk3_20260130.csv \
#     output/combined_LR/combined_LR_all_cousins_n1000_chunk3_20260130.csv \
#     output/unrelated_pool/all_N10000_combined_unrelated_20260702_112423.csv \
#     1 \
#     10 \
#     output/focal_ranking_test/pair_rank_chunks

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  stop(
    "Usage: Rscript code/run_pair_ranking_chunk.R ",
    "<PAIR_FILE> <COMBINED_TRUE_LR_FILE> <UNRELATED_POOL_FILE> ",
    "<CHUNK_INDEX> <CHUNK_SIZE> <OUTPUT_DIR>"
  )
}

PAIR_FILE             <- args[1]
COMBINED_TRUE_LR_FILE <- args[2]
UNRELATED_POOL_FILE   <- args[3]
CHUNK_INDEX           <- as.integer(args[4])
CHUNK_SIZE            <- as.integer(args[5])
OUTPUT_DIR            <- args[6]

setwd("/nfs/turbo/lsa-tlasisi1/podfridge_simulations/")

source("code/module11_ranking_lr_calculator.R")
source("code/module12_ranking_outcome_recorder.R")
source("code/focal_test_helper_fns.R")

kinship_matrix <- fread("data/kinship_coefficients.csv")

tested_relationships <- c("parent_child", "full_siblings")
tested_populations   <- "all"

# -------------------------------------------------------------------------
# Define loci sets from data/core_CODIS_loci.csv
# -------------------------------------------------------------------------

core_loci <- fread("data/core_CODIS_loci.csv")

columns <- c("core_13", "identifiler_15", "expanded_20", "supplementary")

loci_lists <- lapply(columns, function(col) {
  core_loci |>
    filter(get(col) == 1) |>
    pull(locus)
})

names(loci_lists) <- columns

# Add autosomal_29 if needed elsewhere
loci_list <- df_allelefreq |>
  pull(marker) |>
  unique()

loci_lists$autosomal_29 <- loci_list

loci_sets_to_use <- c("core_13", "expanded_20")

loci_sets <- list(
  core_13     = loci_lists$core_13,
  expanded_20 = loci_lists$expanded_20
)

loci_needed <- unique(unlist(loci_sets))

# -------------------------------------------------------------------------
# Output directory
# -------------------------------------------------------------------------

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# -------------------------------------------------------------------------
# Load input files
# -------------------------------------------------------------------------

cat("=============================================================\n")
cat("  Pair ranking chunk job\n")
cat("=============================================================\n")
cat(sprintf("  Pair file             : %s\n", PAIR_FILE))
cat(sprintf("  Combined true LR file : %s\n", COMBINED_TRUE_LR_FILE))
cat(sprintf("  Unrelated pool file   : %s\n", UNRELATED_POOL_FILE))
cat(sprintf("  Chunk index           : %d\n", CHUNK_INDEX))
cat(sprintf("  Chunk size            : %d\n", CHUNK_SIZE))
cat(sprintf("  Output dir            : %s\n", OUTPUT_DIR))
cat("=============================================================\n\n")

pair_data <- fread(
  PAIR_FILE,
  select = c(
    "batch_id",
    "pair_id",
    "population",
    "known_relationship",
    "locus",
    "focal_A1",
    "focal_A2"
  )
)

combined_true_lr_data <- fread(
  COMBINED_TRUE_LR_FILE,
  select = c(
    "batch_id",
    "pair_id",
    "population",
    "known_relationship",
    "loci_set",
    "tested_relationship",
    "tested_population",
    "combined_LR",
    "is_correct_rel",
    "is_correct_pop"
  )
)

unrelated_pool_data <- fread(
  UNRELATED_POOL_FILE,
  select = c(
    "batch_id",
    "database_composition",
    "database_label",
    "individual_id",
    "relationship_to_focal",
    "source_frequency_population",
    "population",
    "locus",
    "A1",
    "A2"
  )
)

# Filter to loci needed for ranking
pair_data <- pair_data |>
  filter(locus %in% loci_needed)

unrelated_pool_data <- unrelated_pool_data |>
  filter(locus %in% loci_needed)

# -------------------------------------------------------------------------
# Determine pair IDs for this chunk
# -------------------------------------------------------------------------

all_pairs <- pair_data |>
  distinct(batch_id, pair_id) |>
  arrange(batch_id, pair_id)

start_idx <- ((CHUNK_INDEX - 1) * CHUNK_SIZE) + 1
end_idx   <- min(CHUNK_INDEX * CHUNK_SIZE, nrow(all_pairs))

if (start_idx > nrow(all_pairs)) {
  stop(sprintf(
    "CHUNK_INDEX %d starts beyond available pairs. There are only %d unique (batch_id, pair_id) pairs.",
    CHUNK_INDEX,
    nrow(all_pairs)
  ))
}

pairs_this_chunk <- all_pairs[start_idx:end_idx, ]

cat(sprintf("Total unique (batch_id, pair_id) pairs in file: %d\n", nrow(all_pairs)))
cat(sprintf("Processing pair index range: %d to %d\n", start_idx, end_idx))
cat(sprintf("Processing %d pairs:\n", nrow(pairs_this_chunk)))
print(pairs_this_chunk)

# -------------------------------------------------------------------------
# Function to rank one pair
# -------------------------------------------------------------------------

rank_one_existing_pair <- function(target_batch_id, target_pair_id) {
  
  cat("\n-------------------------------------------------------------\n")
  cat(sprintf("Processing batch_id: %s, pair_id: %s\n", target_batch_id, target_pair_id))
  cat("-------------------------------------------------------------\n")
  
  pair_data_one_pair <- pair_data |>
    filter(batch_id == target_batch_id, pair_id == target_pair_id)
  
  combined_true_lr_one_pair <- combined_true_lr_data |>
    filter(batch_id == target_batch_id, pair_id == target_pair_id)
  
  if (nrow(pair_data_one_pair) == 0) {
    stop(sprintf("No pair_data rows found for batch_id %s, pair_id %s", target_batch_id, target_pair_id))
  }
  
  if (nrow(combined_true_lr_one_pair) == 0) {
    stop(sprintf("No combined_true_lr rows found for batch_id %s, pair_id %s", target_batch_id, target_pair_id))
  }
  
  pair_meta <- pair_data_one_pair |>
    distinct(batch_id, pair_id, population, known_relationship)
  
  if (nrow(pair_meta) != 1) {
    stop(sprintf(
      "Expected exactly one metadata row for pair_id %s but found %d",
      target_pair_id,
      nrow(pair_meta)
    ))
  }
  
  original_batch_id <- pair_meta$batch_id[1]
  original_pair_id  <- pair_meta$pair_id[1]
  true_relationship <- pair_meta$known_relationship[1]
  true_population   <- pair_meta$population[1]
  
  # Existing true relative combined LR
  true_combined <- format_existing_true_combined_lr(
    combined_lr_one_pair = combined_true_lr_one_pair,
    tested_relationships = tested_relationships,
    tested_populations   = tested_populations,
    loci_sets_to_use     = loci_sets_to_use
  )
  
  # Ensure focal_id exists for Module 12
  true_combined$focal_id <- original_pair_id
  
  # Assemble focal vs unrelated database
  paired_unrelated_db <- assemble_unrelated_database_from_existing_pair(
    pair_data_one_pair  = pair_data_one_pair,
    unrelated_pool_data = unrelated_pool_data,
    batch_id            = paste0("pair_", original_pair_id, "_unrelateds")
  )
  
  # Calculate LRs for focal vs unrelateds
  unrelated_combined <- calculate_ranking_lrs(
    paired_db              = paired_unrelated_db,
    tested_relationships   = tested_relationships,
    tested_populations     = tested_populations,
    loci_sets              = loci_sets,
    allele_frequency_data  = df_allelefreq,
    kinship_coefficients   = kinship_matrix,
    output_dir             = OUTPUT_DIR,
    save_results           = FALSE
  )
  
  unrelated_combined$focal_id <- original_pair_id
  
  # Bind unrelated candidates + true relative
  ranking_input <- bind_rows(
    unrelated_combined,
    true_combined
  )
  
  # Rank true relative
  ranking <- record_ranking_outcomes(
    lr_results   = ranking_input,
    top_n        = 200,
    output_dir   = OUTPUT_DIR,
    save_results = FALSE
  )
  
  outcomes <- ranking$outcomes |>
    mutate(
      original_batch_id = original_batch_id,
      original_pair_id  = original_pair_id,
      true_population   = true_population,
      true_relationship = known_relationship,
      top_200 = rank <= 200,
      top_100 = rank <= 100,
      top_50  = rank <= 50,
      top_10  = rank <= 10
    ) |>
    select(
      original_batch_id,
      original_pair_id,
      true_population,
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
  
  return(outcomes)
}

# -------------------------------------------------------------------------
# Run this chunk
# -------------------------------------------------------------------------

chunk_results_list <- list()
failure_list <- list()

for (i in seq_len(nrow(pairs_this_chunk))) {
  
  bid <- pairs_this_chunk$batch_id[i]
  pid <- pairs_this_chunk$pair_id[i]
  key <- paste(bid, pid, sep = "__")
  
  result <- tryCatch(
    {
      rank_one_existing_pair(bid, pid)
    },
    error = function(e) {
      message(sprintf("ERROR for batch_id %s, pair_id %s: %s", bid, pid, e$message))
      
      failure_list[[key]] <<- data.frame(
        batch_id = bid,
        pair_id = pid,
        error_message = e$message,
        stringsAsFactors = FALSE
      )
      
      NULL
    }
  )
  
  if (!is.null(result)) {
    chunk_results_list[[key]] <- result
  }
  
  # Encourage memory cleanup between pairs
  gc()
}

# -------------------------------------------------------------------------
# Save chunk outputs
# -------------------------------------------------------------------------

safe_pair_file <- tools::file_path_sans_ext(basename(PAIR_FILE))
safe_pair_file <- gsub("[^A-Za-z0-9_\\-]", "_", safe_pair_file)

chunk_out_file <- file.path(
  OUTPUT_DIR,
  sprintf(
    "ranking_outcomes_%s_chunkTask%04d_pairs%05d-%05d.csv",
    safe_pair_file,
    CHUNK_INDEX,
    start_idx,
    end_idx
  )
)

if (length(chunk_results_list) > 0) {
  chunk_results <- rbindlist(chunk_results_list, fill = TRUE)
  fwrite(chunk_results, chunk_out_file)
  
  cat("\nSaved chunk ranking outcomes:\n")
  cat(sprintf("  %s\n", chunk_out_file))
  cat(sprintf("  Rows saved: %d\n", nrow(chunk_results)))
} else {
  warning("No successful pair results in this chunk.")
}

if (length(failure_list) > 0) {
  failures <- rbindlist(failure_list, fill = TRUE)
  
  failure_out_file <- file.path(
    OUTPUT_DIR,
    sprintf(
      "ranking_failures_%s_chunkTask%04d_pairs%05d-%05d.csv",
      safe_pair_file,
      CHUNK_INDEX,
      start_idx,
      end_idx
    )
  )
  
  fwrite(failures, failure_out_file)
  
  cat("\nSaved failure report:\n")
  cat(sprintf("  %s\n", failure_out_file))
}

cat("\nCompleted at:", format(Sys.time()), "\n")
cat("SUCCESS\n")