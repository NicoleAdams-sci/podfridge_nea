##################################
# focal test helper functions    #
##################################

# converts existing combined_LR_...csv row into the format Module 12 expects

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





# creates Module 11 input for unrelated candidates only

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