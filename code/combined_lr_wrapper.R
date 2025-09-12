##### Combined LR Calculation Wrapper Script #####
# USAGE: Rscript code/combined_lr_wrapper.R <LR_CSV_FILE>
# Example: Rscript code/combined_lr_wrapper.R output/LR/LR_AfAm_parent_child_n1000_chunk1_20250905.csv

# Load in arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript combined_lr_wrapper.R <lr_csv_file>")
}

LR_FILE <- args[1]

# Validate input file exists
if (!file.exists(LR_FILE)) {
  stop(paste("Input file does not exist:", LR_FILE))
}

cat("Processing LR file:", LR_FILE, "\n")
cat("Started at:", format(Sys.time()), "\n")

# Load packages
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(tictoc)
})

# Source required modules
source("code/module5_combined_LR.R")

#### Read and validate LR data ####
cat("Reading LR data...\n")
lr_data <- fread(LR_FILE)

# Validate required columns
required_cols <- c("batch_id", "pair_id", "population", "known_relationship", 
                   "locus", "tested_relationship", "tested_population", "LR")
missing_cols <- setdiff(required_cols, names(lr_data))
if (length(missing_cols) > 0) {
  stop(paste("Missing required columns in LR file:", paste(missing_cols, collapse = ", ")))
}

cat("Found", nrow(lr_data), "LR calculations for", 
    length(unique(lr_data$pair_id)), "unique pairs\n")
cat("Loci in data:", length(unique(lr_data$locus)), "loci\n")
cat("Tested relationships:", paste(unique(lr_data$tested_relationship), collapse = ", "), "\n")
cat("Tested populations:", paste(unique(lr_data$tested_population), collapse = ", "), "\n")

#### Calculate combined LRs ####
cat("Calculating combined LRs across loci sets...\n")
tic("Combined LR calculation")

combined_lr_results <- calculate_combined_lr(
  single_locus_results = lr_data,
  loci_sets = loci_lists
)

toc()

cat("Generated", nrow(combined_lr_results), "combined LR calculations\n")
cat("Loci sets:", paste(unique(combined_lr_results$loci_set), collapse = ", "), "\n")

#### Write output immediately ####
cat("Writing output...\n")
tic("Writing output")

# Generate output filename
output_file <- gsub("LR_", "combined_LR_", LR_FILE)
output_file <- gsub("output/LR/", "output/combined_LR/", output_file)

# Create Combined_LR output directory if it doesn't exist
dir.create("output/combined_LR", showWarnings = FALSE, recursive = TRUE)

# Write results immediately
fwrite(combined_lr_results, output_file, quote = FALSE, row.names = FALSE)

# Clear from memory
rm(combined_lr_results, lr_data)

toc()

cat("Output written to:", output_file, "\n")
cat("Completed at:", format(Sys.time()), "\n")
cat("SUCCESS: Combined LR calculation completed for", LR_FILE, "\n")