##### Generate Unrelated Pool #####
# USAGE: Rscript code/generate_unrelated_pool.R <N_UNRELATED_PER_POP>
# Example: Rscript code/generate_unrelated_pool.R 250

# ------------------------------------------------------------------------------
# 0. Parse arguments
# ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript generate_unrelated_pool.R <N_UNRELATED_PER_POP>")
}

N_UNRELATED_PER_POP <- as.integer(args[1])
POPULATIONS         <- c("AfAm", "Cauc", "Hispanic", "Asian")
OUTPUT_DIR          <- "output/unrelated_pool"

cat("Started at:", format(Sys.time()), "\n\n")

# ------------------------------------------------------------------------------
# 1. Load packages and source modules
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(purrr)
})

source("code/LR_kinship_utility_functions.R")
source("code/module1_allele_simulator.R")
source("code/module2_STR_profile_simulator.R")
source("code/module8_unrelated_pool_generator.R")

# ------------------------------------------------------------------------------
# 2. Generate unrelated pool
# ------------------------------------------------------------------------------

cat("Generating unrelated pool...\n")
cat(sprintf("  Populations          : %s\n", paste(POPULATIONS, collapse = ", ")))
cat(sprintf("  Per population       : %d\n", N_UNRELATED_PER_POP))
cat(sprintf("  Total individuals    : %d\n", N_UNRELATED_PER_POP * length(POPULATIONS)))
cat(sprintf("  Output directory     : %s\n\n", OUTPUT_DIR))

results <- generate_multiple_pop_unrelated(
  populations           = POPULATIONS,
  n_unrelated_per_pop   = N_UNRELATED_PER_POP,
  loci_list             = loci_list,
  allele_frequency_data = df_allelefreq,
  output_dir            = OUTPUT_DIR,
  use_single_datetime   = TRUE
)

# ------------------------------------------------------------------------------
# 3. Summary
# ------------------------------------------------------------------------------

cat("\n=============================================================\n")
cat("  Unrelated pool generation complete\n")
cat("=============================================================\n")
cat(sprintf("  Files saved to: %s\n\n", OUTPUT_DIR))
print(results[, c("population", "n_unrelated", "filename")])
cat(sprintf("\n  Total individuals: %d\n", sum(results$n_unrelated)))
cat("Completed at:", format(Sys.time()), "\n")
cat("SUCCESS\n")
