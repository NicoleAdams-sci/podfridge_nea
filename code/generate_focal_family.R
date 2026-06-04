##### Generate Focal Families #####
# USAGE: Rscript code/generate_focal_family.R <N_FOCAL> <POPULATION> <FAMILY_TYPE> [CUSTOM_COUNTS]
#
# Examples:
#   Rscript code/generate_focal_family.R 10 Asian minimal
#   Rscript code/generate_focal_family.R 10 Asian custom "second_cousins=1"
#   Rscript code/generate_focal_family.R 10 Asian custom "parent_child=1,cousins=1,second_cousins=1"

# ------------------------------------------------------------------------------
# 0. Parse arguments
# ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript generate_focal_family.R <N_FOCAL> <POPULATION> <FAMILY_TYPE> [CUSTOM_COUNTS]")
}

N_FOCAL      <- as.integer(args[1])
POPULATION   <- args[2]
FAMILY_TYPE  <- args[3]
CUSTOM_COUNTS <- if (length(args) >= 4) args[4] else ""
OUTPUT_DIR   <- "output/focal_ranking_test"

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
source("code/module3_related_individual_simulator.R")
source("code/module7_single_pop_focal_generator.R")

# Note: must be named kinship_matrix — module 3 references this as a global variable
kinship_matrix <- fread("data/kinship_coefficients.csv")

# ------------------------------------------------------------------------------
# 2. Define family structure
# ------------------------------------------------------------------------------

if (FAMILY_TYPE == "custom") {

  # Parse "rel1=n1,rel2=n2,..." into a named list
  if (CUSTOM_COUNTS == "" || is.na(CUSTOM_COUNTS)) {
    stop("FAMILY_TYPE is 'custom' but no CUSTOM_COUNTS provided.\n",
         "Example: \"second_cousins=1\" or \"parent_child=1,cousins=2\"")
  }

  pairs <- strsplit(CUSTOM_COUNTS, ",")[[1]]
  relationship_counts <- lapply(pairs, function(p) {
    parts <- strsplit(trimws(p), "=")[[1]]
    if (length(parts) != 2)
      stop(sprintf("Invalid custom count format: '%s'. Expected 'relationship=count'", p))
    as.integer(trimws(parts[2]))
  })
  names(relationship_counts) <- sapply(pairs, function(p) trimws(strsplit(p, "=")[[1]][1]))

  # Validate relationship types
  valid_rels <- kinship_matrix$relationship_type
  invalid <- setdiff(names(relationship_counts), valid_rels)
  if (length(invalid) > 0)
    stop(sprintf("Invalid relationship type(s): %s\nValid types: %s",
                 paste(invalid, collapse = ", "),
                 paste(valid_rels, collapse = ", ")))

} else {
  relationship_counts <- create_standard_family_structure(FAMILY_TYPE)
}

cat("Generating focal families...\n")
cat(sprintf("  Population       : %s\n", POPULATION))
cat(sprintf("  N focal          : %d\n", N_FOCAL))
cat(sprintf("  Family type      : %s\n", FAMILY_TYPE))
cat(sprintf("  Relationship structure:\n"))
for (rel in names(relationship_counts)) {
  cat(sprintf("    - %d %s\n", relationship_counts[[rel]], rel))
}
cat(sprintf("  Output directory : %s\n\n", OUTPUT_DIR))

# ------------------------------------------------------------------------------
# 3. Generate focal families
# ------------------------------------------------------------------------------

result <- generate_single_pop_focal(
  population            = POPULATION,
  n_focal               = N_FOCAL,
  relationship_counts   = relationship_counts,
  loci_list             = loci_list,
  allele_frequency_data = df_allelefreq,
  kinship_coefficients  = kinship_matrix,
  output_dir            = OUTPUT_DIR
)

# ------------------------------------------------------------------------------
# 4. Summary
# ------------------------------------------------------------------------------

cat("\n=============================================================\n")
cat("  Focal family generation complete\n")
cat("=============================================================\n")
cat(sprintf("  File saved to: %s\n\n", result$file_path))
cat(sprintf("  Focal individuals : %d\n", N_FOCAL))
cat(sprintf("  Relatives per focal:\n"))
for (rel in names(relationship_counts)) {
  cat(sprintf("    - %d %s\n", relationship_counts[[rel]], rel))
}
cat(sprintf("  Total rows in file: %d\n", nrow(result$data)))
cat("\nCompleted at:", format(Sys.time()), "\n")
cat("SUCCESS\n")
