##### Combine Ranking LR Results #####
# USAGE: Rscript code/combine_ranking_lrs.R <RELATIVE_TYPE>
# Example: Rscript code/combine_ranking_lrs.R parent_child
# Example: Rscript code/combine_ranking_lrs.R second_cousins
#
# Run after all test_module11.sh array tasks complete.
# Reads all lr_focal_<id>_<RELATIVE_TYPE>.csv files from
# output/focal_ranking_test/, stacks them into one combined file.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript combine_ranking_lrs.R <RELATIVE_TYPE>")
}

RELATIVE_TYPE <- args[1]
OUTPUT_DIR    <- "output/focal_ranking_test"

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

cat("=== Combining Ranking LR Results ===\n")
cat(sprintf("Relative type: %s\n", RELATIVE_TYPE))
cat("Started at:", format(Sys.time()), "\n\n")

# Find per-replicate files for this relative type
pattern  <- paste0("^lr_focal_[0-9]+_", RELATIVE_TYPE, "\\.csv$")
lr_files <- list.files(OUTPUT_DIR, pattern = pattern, full.names = TRUE)

if (length(lr_files) == 0)
  stop(sprintf("No lr_focal_*_%s.csv files found in %s", RELATIVE_TYPE, OUTPUT_DIR))

cat(sprintf("Found %d replicate files:\n", length(lr_files)))
for (f in lr_files) cat(sprintf("  %s\n", basename(f)))

# Load, cast combined_LR to numeric, and combine
all_results <- rbindlist(lapply(lr_files, fread))
all_results[, combined_LR := as.numeric(combined_LR)]

cat(sprintf("\nCombined: %s total rows across %d replicates\n",
            format(nrow(all_results), big.mark = ","),
            length(unique(all_results$focal_id))))

# Save combined file
datetime_str <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_file <- file.path(OUTPUT_DIR,
                      sprintf("ranking_lrs_all_replicates_%s_%s.csv",
                              RELATIVE_TYPE, datetime_str))
fwrite(all_results, out_file)
cat(sprintf("Saved: %s\n", out_file))

cat("\nCompleted at:", format(Sys.time()), "\n")
cat("SUCCESS\n")
