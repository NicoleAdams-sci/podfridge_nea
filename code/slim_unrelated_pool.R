##### Slim an Unrelated Pool for the Focal Ranking Test #####
#
# Reads a large *_combined_unrelated_<dt>.csv pool ONCE and writes a small,
# fast-loading .fst copy that keeps only:
#   - the loci in the requested loci set(s), and
#   - the columns the ranking LR actually uses
#     (individual_id, relationship_to_focal, population, locus, A1, A2).
#
# This is a pure filter + reformat: every kept value is carried through
# unchanged, so ranking results are identical to reading the full CSV. The
# only purpose is to stop every ranking array task from re-parsing a multi-GB
# CSV and discarding most of it.
#
# It is deliberately SEPARATE from pool generation: one pool can be slimmed
# for different loci sets, and changing the loci set only means re-running this
# (seconds), never regenerating the pool.
#
# USAGE:
#   Rscript code/slim_unrelated_pool.R <POOL_FILE> <LOCI_SETS> [OUTPUT_DIR]
#
#   <POOL_FILE>   Path to the pool CSV, e.g.
#                 output/unrelated_pool/all_N1000000_combined_unrelated_20260813_132702.csv
#   <LOCI_SETS>   Comma-separated loci-set names from data/core_CODIS_loci.csv
#                 column headers (e.g. core_13,expanded_20). Use these SAME
#                 names in run_focal_ranking.R's loci_sets_to_use.
#   [OUTPUT_DIR]  Optional. Defaults to the pool file's own directory.
#
# EXAMPLES:
#   Rscript code/slim_unrelated_pool.R \
#     output/unrelated_pool/all_N1000000_combined_unrelated_20260813_132702.csv \
#     core_13,expanded_20
#
# The output path is chosen by derive_slim_pool_path() in
# focal_test_helper_fns.R, e.g.
#   output/unrelated_pool/all_N1000000_slim_core13-expanded20.fst
# run_focal_ranking.R reconstructs that same path to find the slim.

suppressPackageStartupMessages({
  library(data.table)
  library(fst)
})

setwd("/nfs/turbo/lsa-tlasisi1/podfridge_simulations/")

source("code/focal_test_helper_fns.R")   # derive_slim_pool_path(), loci_set_tag()

# The columns the ranking path needs from the pool (see
# assemble_unrelated_database_from_existing_pair()). Everything else
# (batch_id, database_composition, database_label, source_frequency_population)
# is unused by ranking and dropped.
POOL_COLS <- c("individual_id", "relationship_to_focal", "population",
               "locus", "A1", "A2")

# ---------------------------------------------------------------------------
# 1. Parse arguments
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    "Usage: Rscript code/slim_unrelated_pool.R <POOL_FILE> <LOCI_SETS> [OUTPUT_DIR]\n",
    "  e.g. Rscript code/slim_unrelated_pool.R ",
    "output/unrelated_pool/all_N1000000_combined_unrelated_<dt>.csv core_13,expanded_20",
    call. = FALSE
  )
}

POOL_FILE  <- args[1]
LOCI_SETS  <- trimws(strsplit(args[2], ",", fixed = TRUE)[[1]])
OUTPUT_DIR <- if (length(args) >= 3) args[3] else NULL

if (!file.exists(POOL_FILE)) {
  stop("Pool file does not exist: ", POOL_FILE, call. = FALSE)
}

# ---------------------------------------------------------------------------
# 2. Derive loci_needed the SAME way run_focal_ranking.R does
#    (union of the requested loci sets from data/core_CODIS_loci.csv).
# ---------------------------------------------------------------------------

core_loci <- fread("data/core_CODIS_loci.csv")

resolve_loci_set <- function(set_name) {
  if (set_name == "autosomal_29") {
    # "all loci" — no locus filter; keep whatever is in the pool.
    return(NA_character_)  # sentinel handled below
  }
  if (!set_name %in% names(core_loci)) {
    stop("Unknown loci set '", set_name, "'. Available: ",
         paste(setdiff(names(core_loci), "locus"), collapse = ", "),
         " (or autosomal_29 to keep all loci)", call. = FALSE)
  }
  core_loci$locus[core_loci[[set_name]] == 1]
}

resolved <- lapply(LOCI_SETS, resolve_loci_set)
keep_all_loci <- any(vapply(resolved, function(x) any(is.na(x)), logical(1)))
loci_needed <- if (keep_all_loci) NULL else sort(unique(unlist(resolved)))

# ---------------------------------------------------------------------------
# 3. Read the pool once (only the needed columns, forced to character so
#    allele/key values round-trip exactly as written by Module 8).
# ---------------------------------------------------------------------------

cat("=============================================================\n")
cat("  Slim unrelated pool\n")
cat("=============================================================\n")
cat(sprintf("  Pool file   : %s\n", POOL_FILE))
cat(sprintf("  Loci sets   : %s\n", paste(LOCI_SETS, collapse = ", ")))
cat(sprintf("  Loci kept   : %s\n",
            if (keep_all_loci) "ALL (autosomal_29)" else paste(loci_needed, collapse = ", ")))
cat(sprintf("  Started at  : %s\n", format(Sys.time())))
cat("=============================================================\n\n")

pool <- fread(
  POOL_FILE,
  select     = POOL_COLS,
  colClasses = setNames(rep("character", length(POOL_COLS)), POOL_COLS)
)

n_rows_in <- nrow(pool)

# ---------------------------------------------------------------------------
# 4. Filter to loci_needed (unless keeping all) and sanity-check
# ---------------------------------------------------------------------------

if (!keep_all_loci) {
  present <- unique(pool$locus)
  missing_loci <- setdiff(loci_needed, present)
  if (length(missing_loci) > 0) {
    stop("Pool is missing requested loci: ", paste(missing_loci, collapse = ", "),
         "\n(Pool contains: ", paste(sort(present), collapse = ", "), ")",
         call. = FALSE)
  }
  pool <- pool[locus %in% loci_needed]
}

n_rows_out <- nrow(pool)
n_ind      <- length(unique(pool$individual_id))
n_loci     <- length(unique(pool$locus))

# ---------------------------------------------------------------------------
# 5. Write the slim .fst
# ---------------------------------------------------------------------------

slim_path <- derive_slim_pool_path(POOL_FILE, LOCI_SETS, output_dir = OUTPUT_DIR)
dir.create(dirname(slim_path), recursive = TRUE, showWarnings = FALSE)

write_fst(pool, slim_path, compress = 60)

size_in_gb  <- file.size(POOL_FILE) / 1024^3
size_out_gb <- file.size(slim_path) / 1024^3

cat("\n=============================================================\n")
cat("  Slim complete\n")
cat("=============================================================\n")
cat(sprintf("  Slim file        : %s\n", slim_path))
cat(sprintf("  Individuals      : %d\n", n_ind))
cat(sprintf("  Loci             : %d\n", n_loci))
cat(sprintf("  Rows in  -> out  : %d -> %d\n", n_rows_in, n_rows_out))
cat(sprintf("  Size  in -> out  : %.2f GB -> %.2f GB\n", size_in_gb, size_out_gb))
cat(sprintf("  Completed at     : %s\n", format(Sys.time())))
cat("SUCCESS\n")
