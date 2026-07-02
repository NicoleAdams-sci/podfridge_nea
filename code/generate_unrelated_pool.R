##### Generate Unrelated Pool / Database #####
# 20260625
#
# USAGE:
#   Rscript code/generate_unrelated_pool.R all <N_TOTAL>
#   Rscript code/generate_unrelated_pool.R single <N_TOTAL> <POPULATION>
#   Rscript code/generate_unrelated_pool.R mixed-proportions <N_TOTAL> <POP1=W1,POP2=W2,...>
#   Rscript code/generate_unrelated_pool.R mixed-counts <POP1=N1,POP2=N2,...>
#
# EXAMPLES:
#   Rscript code/generate_unrelated_pool.R all 100000
#   Rscript code/generate_unrelated_pool.R single 100000 Cauc
#   Rscript code/generate_unrelated_pool.R mixed-proportions 100000 AfAm=0.15,Cauc=0.55,Hispanic=0.20,Asian=0.10
#   Rscript code/generate_unrelated_pool.R mixed-counts AfAm=15000,Cauc=55000,Hispanic=20000,Asian=10000

# ------------------------------------------------------------------------------
# 0. Helper functions
# ------------------------------------------------------------------------------

print_usage_and_stop <- function() {
  stop(
    paste0(
      "Usage:\n",
      "  Rscript code/generate_unrelated_pool.R all <N_TOTAL>\n",
      "  Rscript code/generate_unrelated_pool.R single <N_TOTAL> <POPULATION>\n",
      "  Rscript code/generate_unrelated_pool.R mixed-proportions <N_TOTAL> <POP1=W1,POP2=W2,...>\n",
      "  Rscript code/generate_unrelated_pool.R mixed-counts <POP1=N1,POP2=N2,...>\n",
      "\nExamples:\n",
      "  Rscript code/generate_unrelated_pool.R all 100000\n",
      "  Rscript code/generate_unrelated_pool.R single 100000 Cauc\n",
      "  Rscript code/generate_unrelated_pool.R mixed-proportions 100000 AfAm=0.15,Cauc=0.55,Hispanic=0.20,Asian=0.10\n",
      "  Rscript code/generate_unrelated_pool.R mixed-counts AfAm=15000,Cauc=55000,Hispanic=20000,Asian=10000\n"
    ),
    call. = FALSE
  )
}

parse_named_numeric <- function(x) {
  # Parses strings like:
  #   "AfAm=0.15,Cauc=0.55,Hispanic=0.20,Asian=0.10"
  # or
  #   "AfAm=15000,Cauc=55000,Hispanic=20000,Asian=10000"
  
  pieces <- strsplit(x, ",", fixed = TRUE)[[1]]
  pieces <- trimws(pieces)
  
  values <- numeric(length(pieces))
  names_out <- character(length(pieces))
  
  for (i in seq_along(pieces)) {
    kv <- strsplit(pieces[i], "=", fixed = TRUE)[[1]]
    
    if (length(kv) != 2) {
      stop(
        "Invalid named numeric specification: ", pieces[i],
        "\nExpected format like AfAm=0.15,Cauc=0.55 or AfAm=15000,Cauc=55000",
        call. = FALSE
      )
    }
    
    names_out[i] <- trimws(kv[1])
    values[i] <- as.numeric(trimws(kv[2]))
  }
  
  if (any(is.na(values))) {
    stop("Could not parse one or more numeric values from: ", x, call. = FALSE)
  }
  
  names(values) <- names_out
  values
}

allocate_counts_from_proportions <- function(proportions, n_total) {
  # Converts proportions/weights to exact integer counts summing to n_total.
  # Uses largest remainder method.
  
  if (is.null(names(proportions)) || any(names(proportions) == "")) {
    stop("Proportions must be named, e.g. AfAm=0.15,Cauc=0.55", call. = FALSE)
  }
  
  if (any(is.na(proportions)) || any(proportions < 0)) {
    stop("All proportions/weights must be non-missing and >= 0", call. = FALSE)
  }
  
  if (sum(proportions) <= 0) {
    stop("At least one proportion/weight must be > 0", call. = FALSE)
  }
  
  if (is.na(n_total) || n_total < 1) {
    stop("N_TOTAL must be a positive integer", call. = FALSE)
  }
  
  # Normalize so user can pass either true proportions summing to 1
  # or weights such as Cauc=55,AfAm=15,Hispanic=20,Asian=10.
  proportions <- proportions / sum(proportions)
  
  raw_counts <- proportions * n_total
  base_counts <- floor(raw_counts)
  
  remainder <- n_total - sum(base_counts)
  
  if (remainder > 0) {
    fractional_parts <- raw_counts - base_counts
    add_to <- order(fractional_parts, decreasing = TRUE)[seq_len(remainder)]
    base_counts[add_to] <- base_counts[add_to] + 1
  }
  
  counts <- as.integer(base_counts)
  names(counts) <- names(proportions)
  
  counts
}

relabel_individual_ids <- function(df, start_index = 1) {
  # Module 8 starts each generated population at unrel_001.
  # For combined/mixed databases, this relabels individuals globally:
  # unrel_000001, unrel_000002, ...
  
  old_ids <- unique(df$individual_id)
  new_ids <- paste0("unrel_", sprintf("%06d", seq(
    from = start_index,
    length.out = length(old_ids)
  )))
  
  id_map <- data.frame(
    individual_id = old_ids,
    new_individual_id = new_ids,
    stringsAsFactors = FALSE
  )
  
  df <- dplyr::left_join(df, id_map, by = "individual_id")
  
  df$individual_id <- df$new_individual_id
  df$new_individual_id <- NULL
  
  df
}

make_safe_label <- function(x) {
  gsub("[^A-Za-z0-9_\\-]", "_", x)
}

# ------------------------------------------------------------------------------
# 1. Parse arguments
# ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  print_usage_and_stop()
}

MODE <- args[1]
OUTPUT_DIR <- "output/unrelated_pool"

valid_modes <- c(
  "all",
  "single",
  "mixed-proportions",
  "mixed-counts"
)

if (!MODE %in% valid_modes) {
  stop(
    "Invalid mode: ", MODE,
    "\nValid modes are: ", paste(valid_modes, collapse = ", "),
    call. = FALSE
  )
}

cat("Started at:", format(Sys.time()), "\n\n")

# ------------------------------------------------------------------------------
# 2. Load packages and source modules
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
# 3. Convert mode into exact counts by population
# ------------------------------------------------------------------------------

if (MODE == "all") {
  
  if (length(args) != 2) {
    print_usage_and_stop()
  }
  
  n_total <- as.integer(args[2])
  
  if (is.na(n_total) || n_total < 1) {
    stop("N_TOTAL must be a positive integer", call. = FALSE)
  }
  
  counts_by_population <- c(all = n_total)
  database_label <- paste0("all_N", n_total)
  
} else if (MODE == "single") {
  
  if (length(args) != 3) {
    print_usage_and_stop()
  }
  
  n_total <- as.integer(args[2])
  pop <- args[3]
  
  if (is.na(n_total) || n_total < 1) {
    stop("N_TOTAL must be a positive integer", call. = FALSE)
  }
  
  counts_by_population <- setNames(n_total, pop)
  database_label <- paste0("single_", pop, "_N", n_total)
  
} else if (MODE == "mixed-proportions") {
  
  if (length(args) != 3) {
    print_usage_and_stop()
  }
  
  n_total <- as.integer(args[2])
  proportions <- parse_named_numeric(args[3])
  
  counts_by_population <- allocate_counts_from_proportions(
    proportions = proportions,
    n_total = n_total
  )
  
  database_label <- paste0(
    "mixed_proportions_N",
    n_total,
    "_",
    paste(names(proportions), proportions, sep = "-", collapse = "_")
  )
  
} else if (MODE == "mixed-counts") {
  
  if (length(args) != 2) {
    print_usage_and_stop()
  }
  
  counts_by_population <- parse_named_numeric(args[2])
  counts_by_population <- as.integer(counts_by_population)
  
  if (any(is.na(counts_by_population)) || any(counts_by_population < 0)) {
    stop("All mixed-counts values must be non-missing integers >= 0", call. = FALSE)
  }
  
  counts_by_population <- counts_by_population[counts_by_population > 0]
  
  if (length(counts_by_population) == 0) {
    stop("At least one population count must be > 0", call. = FALSE)
  }
  
  database_label <- paste0(
    "mixed_counts_",
    paste(names(counts_by_population), counts_by_population, sep = "-", collapse = "_")
  )
}

# ------------------------------------------------------------------------------
# 4. Validate requested populations
# ------------------------------------------------------------------------------

available_populations <- unique(df_allelefreq$population)
missing_populations <- setdiff(names(counts_by_population), available_populations)

if (length(missing_populations) > 0) {
  stop(
    "The following requested populations are not present in df_allelefreq$population: ",
    paste(missing_populations, collapse = ", "),
    "\nAvailable populations are: ",
    paste(available_populations, collapse = ", "),
    call. = FALSE
  )
}

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

shared_datetime <- format(Sys.time(), "%Y%m%d_%H%M%S")

cat("=============================================================\n")
cat("  Generate Unrelated Database\n")
cat("=============================================================\n")
cat(sprintf("  Mode             : %s\n", MODE))
cat(sprintf("  Database label   : %s\n", database_label))
cat(sprintf("  Output directory : %s\n", OUTPUT_DIR))
cat(sprintf("  Datetime         : %s\n", shared_datetime))
cat("  Counts by population:\n")
print(counts_by_population)
cat(sprintf("  Total individuals: %d\n", sum(counts_by_population)))
cat("=============================================================\n\n")

# ------------------------------------------------------------------------------
# 5. Generate each requested population using existing Module 8 function
# ------------------------------------------------------------------------------

all_chunks <- list()
summary_rows <- list()
next_start_index <- 1

for (pop in names(counts_by_population)) {
  
  n_pop <- counts_by_population[[pop]]
  
  cat("\n-------------------------------------------------------------\n")
  cat(sprintf("Generating population chunk: %s, N = %d\n", pop, n_pop))
  cat("-------------------------------------------------------------\n")
  
  result <- generate_unrelated_pool(
    population            = pop,
    n_unrelated           = n_pop,
    loci_list             = loci_list,
    allele_frequency_data = df_allelefreq,
    output_dir            = OUTPUT_DIR,
    custom_datetime       = shared_datetime
  )
  
  chunk <- result$data
  
  # Add metadata useful for downstream tracking.
  chunk$database_composition <- MODE
  chunk$database_label <- database_label
  chunk$source_frequency_population <- pop
  
  # Relabel individual IDs so combined file has globally unique IDs.
  chunk <- relabel_individual_ids(
    df = chunk,
    start_index = next_start_index
  )
  
  n_generated_individuals <- length(unique(chunk$individual_id))
  next_start_index <- next_start_index + n_generated_individuals
  
  all_chunks[[pop]] <- chunk
  
  summary_rows[[pop]] <- data.frame(
    population = pop,
    n_unrelated = n_pop,
    module8_file_path = result$file_path,
    module8_filename = basename(result$file_path),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# 6. Combine chunks and save one final database file
# ------------------------------------------------------------------------------

combined_data <- dplyr::bind_rows(all_chunks)

# Reorder columns, keeping any unexpected extra columns at the end.
preferred_cols <- c(
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

existing_preferred_cols <- preferred_cols[preferred_cols %in% names(combined_data)]
other_cols <- setdiff(names(combined_data), existing_preferred_cols)

combined_data <- combined_data[, c(existing_preferred_cols, other_cols)]

safe_label <- make_safe_label(database_label)
combined_filename <- paste0(safe_label, "_combined_unrelated_", shared_datetime, ".csv")
combined_file_path <- file.path(OUTPUT_DIR, combined_filename)

cat("\nSaving combined unrelated database to:\n")
cat(sprintf("  %s\n", combined_file_path))

data.table::fwrite(combined_data, combined_file_path)

summary_df <- dplyr::bind_rows(summary_rows)

# Also save a small composition summary.
summary_filename <- paste0(safe_label, "_composition_summary_", shared_datetime, ".csv")
summary_file_path <- file.path(OUTPUT_DIR, summary_filename)

summary_out <- data.frame(
  database_label = database_label,
  mode = MODE,
  population = names(counts_by_population),
  n_unrelated = as.integer(counts_by_population),
  proportion = as.numeric(counts_by_population) / sum(counts_by_population),
  stringsAsFactors = FALSE
)

data.table::fwrite(summary_out, summary_file_path)

# ------------------------------------------------------------------------------
# 7. Final summary
# ------------------------------------------------------------------------------

cat("\n=============================================================\n")
cat("  Unrelated database generation complete\n")
cat("=============================================================\n")
cat(sprintf("  Combined file      : %s\n", combined_file_path))
cat(sprintf("  Summary file       : %s\n", summary_file_path))
cat(sprintf("  Total individuals  : %d\n", sum(counts_by_population)))
cat(sprintf("  Total rows         : %d\n", nrow(combined_data)))
cat("\n  Composition:\n")
print(summary_out[, c("population", "n_unrelated", "proportion")])
cat("\n  Module 8 intermediate files:\n")
print(summary_df[, c("population", "n_unrelated", "module8_filename")])
cat("\nCompleted at:", format(Sys.time()), "\n")
cat("SUCCESS\n")