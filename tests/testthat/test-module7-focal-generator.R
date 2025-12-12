test_that("generate_single_pop_focal produces expected structure and file output", {
  tmpdir <- tempfile("focal_single")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)
  
  with_mocked_focal_functions({
    rel_counts <- list(parent_child = 1, full_siblings = 2)
    result <- generate_single_pop_focal(
      population = "AfAm",
      n_focal = 2,
      relationship_counts = rel_counts,
      loci_list = mock_loci,
      allele_frequency_data = mock_allele_freq,
      kinship_coefficients = mock_kinship_coeffs,
      output_dir = tmpdir,
      custom_datetime = "20250101_000000"
    )
    
    expect_true(all(expected_focal_columns %in% names(result$data)))
    expect_equal(length(unique(result$data$focal_id)), 2)
    
    individuals_per_family <- 1 + sum(unlist(rel_counts))
    expected_rows <- 2 * individuals_per_family * length(mock_loci)
    expect_equal(nrow(result$data), expected_rows)
    
    expect_setequal(unique(result$data$relationship_to_focal), c("self", names(rel_counts)))
    expect_true(file.exists(result$file_path))
    expect_match(basename(result$file_path), "AfAm_focal_parent_child1_full_siblings2_20250101_000000.csv")
  })
})

test_that("generate_single_pop_focal validates populations and relationship counts", {
  with_mocked_focal_functions({
    expect_error(
      generate_single_pop_focal(
        population = "Martian",
        n_focal = 1,
        relationship_counts = list(parent_child = 1),
        loci_list = mock_loci,
        allele_frequency_data = mock_allele_freq,
        kinship_coefficients = mock_kinship_coeffs
      ),
      "Population .* not found",
      fixed = FALSE
    )
    
    expect_error(
      generate_single_pop_focal(
        population = "AfAm",
        n_focal = 1,
        relationship_counts = list(grandparent = 1),
        loci_list = mock_loci,
        allele_frequency_data = mock_allele_freq,
        kinship_coefficients = mock_kinship_coeffs
      ),
      "Invalid relationship types",
      fixed = TRUE
    )
    
    expect_error(
      generate_single_pop_focal(
        population = "AfAm",
        n_focal = 1,
        relationship_counts = list(parent_child = -1),
        loci_list = mock_loci,
        allele_frequency_data = mock_allele_freq,
        kinship_coefficients = mock_kinship_coeffs
      ),
      "non-negative integers",
      fixed = FALSE
    )
  })
})

test_that("generate_families_by_relationships assigns metadata consistently", {
  with_mocked_focal_functions({
    families <- generate_families_by_relationships(
      n_focal = 1,
      population = "AfAm",
      relationship_counts = list(parent_child = 1, half_siblings = 1, cousins = 0),
      loci_list = mock_loci,
      allele_frequency_data = mock_allele_freq,
      kinship_coefficients = mock_kinship_coeffs,
      batch_id = "batch123"
    )
    
    expect_equal(unique(families$batch_id), "batch123")
    expect_equal(unique(families$family_id), "fam_001")
    expect_equal(unique(families$focal_id), "001")
    expect_setequal(unique(families$relationship_to_focal), c("self", "parent_child", "half_siblings"))
    
    parent_ids <- unique(families$individual_id[families$relationship_to_focal == "parent_child"])
    expect_true(all(grepl("^parentchild1_001$", parent_ids)))
    
    half_sib_ids <- unique(families$individual_id[families$relationship_to_focal == "half_siblings"])
    expect_true(all(grepl("^halfsiblings1_001$", half_sib_ids)))
  })
})

test_that("generate_multiple_pop_focal returns per-pop summaries with shared datetime", {
  tmpdir <- tempfile("focal_multi")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)
  
  with_mocked_focal_functions({
    structures <- list(
      AfAm = list(parent_child = 1),
      Cauc = list(full_siblings = 1)
    )
    
    results <- generate_multiple_pop_focal(
      population_structures = structures,
      n_focal_per_pop = 1,
      loci_list = mock_loci,
      allele_frequency_data = mock_allele_freq,
      kinship_coefficients = mock_kinship_coeffs,
      output_dir = tmpdir,
      use_single_datetime = TRUE
    )
    
    expect_equal(nrow(results), length(structures))
    expect_equal(length(unique(results$datetime_used)), 1)
    expect_true(all(file.exists(results$file_path)))
    expect_setequal(results$population, names(structures))
  })
})
