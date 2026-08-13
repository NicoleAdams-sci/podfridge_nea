##################################
# focal test helper functions    #
##################################

# converts existing combined_LR_...csv row into the format Module 11 expects

format_existing_true_combined_lr <- function(combined_lr_one_pair,
                                             tested_relationships = c("parent_child", "full_siblings"),
                                             tested_populations = "all",
                                             loci_sets_to_use = c("core_13", "core_20")) {
  
  required_cols <- c(
    "batch_id",
    "pair_id",
    "population",
    "known_relationship",
    "loci_set",
    "tested_relationship",
    "tested_population",
    "combined_LR"
  )
  
  missing_cols <- setdiff(required_cols, names(combined_lr_one_pair))
  if (length(missing_cols) > 0) {
    stop("combined_lr_one_pair missing columns: ", paste(missing_cols, collapse = ", "))
  }
  
  pair_ids <- unique(combined_lr_one_pair$pair_id)
  if (length(pair_ids) != 1) {
    stop("combined_lr_one_pair must contain exactly one pair_id")
  }
  
  pair_id_value <- pair_ids[1]
  
  true_rows <- combined_lr_one_pair |>
    dplyr::filter(
      tested_relationship %in% tested_relationships,
      tested_population %in% tested_populations,
      loci_set %in% loci_sets_to_use
    ) |>
    dplyr::mutate(
      combined_LR = as.numeric(combined_LR),
      focal_id = pair_id_value,
      individual_id = paste0("true_relative_", pair_id_value),
      is_true_relative = TRUE
    ) |>
    dplyr::select(
      focal_id,
      pair_id,
      individual_id,
      loci_set,
      tested_relationship,
      tested_population,
      known_relationship,
      combined_LR,
      is_true_relative,
      dplyr::everything()
    )
  
  if (nrow(true_rows) == 0) {
    stop("No true-relative combined LR rows remain after filtering.")
  }
  
  return(true_rows)
}





# creates Module 10 input for unrelated candidates only

assemble_unrelated_database_from_existing_pair <- function(pair_data_one_pair,
                                                           unrelated_pool_data,
                                                           batch_id = NULL) {
  
  required_pair_cols <- c(
    "batch_id",
    "pair_id",
    "population",
    "known_relationship",
    "locus",
    "focal_A1",
    "focal_A2"
  )
  
  required_unrel_cols <- c(
    "individual_id",
    "relationship_to_focal",
    "population",
    "locus",
    "A1",
    "A2"
  )
  
  missing_pair <- setdiff(required_pair_cols, names(pair_data_one_pair))
  missing_unrel <- setdiff(required_unrel_cols, names(unrelated_pool_data))
  
  if (length(missing_pair) > 0) {
    stop("pair_data_one_pair missing columns: ", paste(missing_pair, collapse = ", "))
  }
  
  if (length(missing_unrel) > 0) {
    stop("unrelated_pool_data missing columns: ", paste(missing_unrel, collapse = ", "))
  }
  
  original_pair_id <- unique(pair_data_one_pair$pair_id)
  if (length(original_pair_id) != 1) {
    stop("pair_data_one_pair must contain exactly one pair_id.")
  }
  
  focal_population <- unique(pair_data_one_pair$population)
  if (length(focal_population) != 1) {
    stop("pair_data_one_pair must contain exactly one population.")
  }
  
  if (is.null(batch_id)) {
    batch_id <- paste0("pair_", original_pair_id, "_unrelateds")
  }
  
  focal_genotype <- pair_data_one_pair |>
    dplyr::select(locus, focal_A1, focal_A2) |>
    dplyr::distinct()
  
  unrelated_genotype <- unrelated_pool_data |>
    dplyr::transmute(
      locus,
      ind2_A1 = as.character(A1),
      ind2_A2 = as.character(A2),
      individual_id = paste0("unrelated_", individual_id),
      known_relationship = "unrelated",
      is_true_relative = FALSE
    )
  
  paired_db <- unrelated_genotype |>
    dplyr::left_join(focal_genotype, by = "locus") |>
    dplyr::mutate(
      batch_id = batch_id,
      original_pair_id = original_pair_id,
      population = focal_population
    )
  
  candidate_map <- data.frame(
    individual_id = unique(paired_db$individual_id),
    pair_id = sprintf("unrel_%05d", seq_along(unique(paired_db$individual_id))),
    stringsAsFactors = FALSE
  )
  
  paired_db <- paired_db |>
    dplyr::left_join(candidate_map, by = "individual_id") |>
    dplyr::select(
      batch_id,
      original_pair_id,
      pair_id,
      individual_id,
      population,
      locus,
      focal_A1,
      focal_A2,
      ind2_A1,
      ind2_A2,
      known_relationship,
      is_true_relative
    )
  
  return(paired_db)
}



# ---------------------------------------------------------------------------
# Slim-pool naming convention (single source of truth)
#
# The ranking test can load a pre-slimmed .fst copy of the unrelated pool,
# filtered to just the loci sets in use and the columns the LR actually needs.
# The slim file's name encodes BOTH the pool identity and the loci sets, so
# that (a) two different loci-set slims of the same pool don't overwrite each
# other, and (b) run_focal_ranking.R can locate the correct slim for its
# analysis purely by convention.
#
# slim_unrelated_pool.R WRITES to this path; run_focal_ranking.R READS from it.
# Both call these two functions so the convention lives in exactly one place.
# ---------------------------------------------------------------------------

#' Order-independent, filesystem-safe tag for a set of loci-set names.
#' e.g. c("expanded_20", "core_13") -> "core13-expanded20"
loci_set_tag <- function(loci_set_names) {
  cleaned <- gsub("[^A-Za-z0-9]", "", loci_set_names)
  paste(sort(unique(cleaned)), collapse = "-")
}

#' Derive the slim-pool .fst path for a given pool file and loci sets.
#'
#' @param pool_file       Path to the *_combined_unrelated_<dt>.csv pool file.
#' @param loci_set_names  Character vector of loci-set names (e.g.
#'                        c("core_13", "expanded_20")) — the SAME names used to
#'                        build the slim and the same ones run_focal_ranking.R
#'                        ranks on.
#' @param output_dir      Directory for the slim file. Defaults to the pool
#'                        file's own directory (slim sits next to the pool).
#' @return Character path like
#'         output/unrelated_pool/all_N1000000_slim_core13-expanded20.fst
derive_slim_pool_path <- function(pool_file, loci_set_names, output_dir = NULL) {
  pool_base <- basename(pool_file)

  # Recover the database label (e.g. "all_N1000000") by stripping the
  # "_combined_unrelated_<YYYYMMDD>_<HHMMSS>.csv" (or .fst) suffix.
  database_label <- sub("_combined_unrelated_[0-9]{8}_[0-9]{6}\\.(csv|fst)$", "",
                        pool_base)
  if (identical(database_label, pool_base)) {
    # Name didn't match the expected pattern; fall back to sans-extension.
    database_label <- tools::file_path_sans_ext(pool_base)
  }

  slim_name <- paste0(database_label, "_slim_", loci_set_tag(loci_set_names), ".fst")

  if (is.null(output_dir)) output_dir <- dirname(pool_file)
  file.path(output_dir, slim_name)
}
