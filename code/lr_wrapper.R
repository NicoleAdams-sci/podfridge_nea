##### LR Calculation Wrapper Script #####
# USAGE: Rscript code/lr_wrapper.R <PAIRS_CSV_FILE>
# Example: Rscript code/lr_wrapper.R output/pairs_AfAm_parent_child_n1000_chunk1_20250905.csv

# Load in arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript lr_wrapper.R <pairs_csv_file>")
}

PAIRS_FILE <- args[1]

# Validate input file exists
if (!file.exists(PAIRS_FILE)) {
  stop(paste("Input file does not exist:", PAIRS_FILE))
}

cat("Processing pairs file:", PAIRS_FILE, "\n")
cat("Started at:", format(Sys.time()), "\n")

# Load packages
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(tictoc)
})

# Source required modules
source("code/LR_kinship_utility_functions.R")
source("code/module4_single_locus_LR.R")

#### Load reference datasets (load once, use for all calculations) ####
cat("Loading reference data...\n")

# Load kinship coefficients
kinship_matrix <- fread("data/kinship_coefficients.csv")

# Define testing parameters (fixed for all calculations)
tested_relationships <- c("parent_child", "full_siblings")
tested_populations <- c("AfAm", "Cauc", "Hispanic", "Asian")

#### Read and validate pairs data ####
cat("Reading pairs data...\n")
pairs_data <- fread(PAIRS_FILE)

# Validate required columns
required_cols <- c("batch_id", "pair_id", "population", "locus", "focal_A1", "focal_A2", 
                   "ind2_A1", "ind2_A2", "known_relationship")
missing_cols <- setdiff(required_cols, names(pairs_data))
if (length(missing_cols) > 0) {
  stop(paste("Missing required columns in pairs file:", paste(missing_cols, collapse = ", ")))
}

cat("Found", nrow(pairs_data), "rows with", length(unique(pairs_data$pair_id)), "unique pairs\n")

#### Calculate LRs for each tested relationship ####
all_lr_results <- list()

for (test_rel in tested_relationships) {
  cat("Calculating LRs for tested relationship:", test_rel, "\n")
  tic(paste("LR calculation for", test_rel))
  
  # Calculate LRs for this relationship
  lr_results <- calculate_single_locus_lr(
    pair_data = pairs_data,
    tested_relationship = test_rel,
    tested_populations = tested_populations,
    allele_frequency_data = df_allelefreq,
    kinship_coefficients = kinship_matrix
  )
  
  # Store results
  all_lr_results[[test_rel]] <- lr_results
  
  toc()
  cat("  Generated", nrow(lr_results), "LR calculations\n")
}

#### Combine all results and write immediately ####
cat("Combining results and writing output...\n")
tic("Writing output")

# Combine all LR results
final_lr_results <- bind_rows(all_lr_results)

# Generate output filename
output_file <- gsub("pairs_", "LR_", PAIRS_FILE)
output_file <- gsub("output/", "output/LR/", output_file)

# Create LR output directory if it doesn't exist
dir.create("output/LR", showWarnings = FALSE, recursive = TRUE)

# Write results immediately
fwrite(final_lr_results, output_file, quote = FALSE, row.names = FALSE)

# Clear from memory
rm(final_lr_results, all_lr_results)

toc()

cat("Output written to:", output_file, "\n")
cat("Completed at:", format(Sys.time()), "\n")
cat("SUCCESS: LR calculation completed for", PAIRS_FILE, "\n")