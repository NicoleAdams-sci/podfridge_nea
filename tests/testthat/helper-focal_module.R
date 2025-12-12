# Helper utilities for testing Module 7 focal generators

# Source the module with chdir=TRUE so its internal relative sources resolve
source_module_file <- function(relative_path) {
  candidates <- c(
    relative_path,
    file.path("..", relative_path),
    file.path("..", "..", relative_path)
  )
  target <- candidates[file.exists(candidates)][1]
  if (is.na(target)) {
    stop("Unable to locate ", relative_path, call. = FALSE)
  }
  source(target, chdir = TRUE)
}

if (!exists("generate_single_pop_focal")) {
  source_module_file("code/module7_single_pop_focal_generator.R")
}

mock_loci <- c("D1S1358", "D2S441")

mock_allele_freq <- data.table::data.table(
  population = rep(c("AfAm", "Cauc"), each = 4),
  marker = rep(mock_loci, times = 4),
  allele = rep(c("10", "11"), times = 4),
  frequency = 0.5
)

mock_kinship_coeffs <- data.table::data.table(
  relationship_type = c("parent_child", "full_siblings", "half_siblings", "cousins", "second_cousins", "unrelated"),
  k0 = c(0, 0.25, 0.5, 0.875, 0.9375, 1),
  k1 = c(1, 0.5, 0.5, 0.125, 0.0625, 0),
  k2 = c(0, 0.25, 0, 0, 0, 0)
)

expected_focal_columns <- c("batch_id", "family_id", "focal_id", "individual_id",
                            "relationship_to_focal", "population", "locus", "A1", "A2")

with_mocked_focal_functions <- function(code) {
  env <- globalenv()
  if (!exists("simulate_str_profile", envir = env) ||
      !exists("simulate_related_individual", envir = env)) {
    stop("Simulation functions must exist before mocking.", call. = FALSE)
  }
  original_profile <- get("simulate_str_profile", envir = env)
  original_related <- get("simulate_related_individual", envir = env)
  
  mock_profile <- function(loci_list, population, allele_frequency_data) {
    data.frame(
      population = population,
      locus = loci_list,
      A1 = paste0(loci_list, "_A1"),
      A2 = paste0(loci_list, "_A2"),
      stringsAsFactors = FALSE
    )
  }
  
  mock_related <- function(focal_profile, known_relationship, allele_frequency_data, individual_id) {
    data.frame(
      individual_id = individual_id,
      population = focal_profile$population,
      locus = focal_profile$locus,
      A1 = paste0(individual_id, "_", focal_profile$locus, "_A1"),
      A2 = paste0(individual_id, "_", focal_profile$locus, "_A2"),
      relationship_to_focal = known_relationship,
      stringsAsFactors = FALSE
    )
  }
  
  assign("simulate_str_profile", mock_profile, envir = env)
  assign("simulate_related_individual", mock_related, envir = env)
  
  on.exit({
    assign("simulate_str_profile", original_profile, envir = env)
    assign("simulate_related_individual", original_related, envir = env)
  }, add = TRUE)
  
  eval(substitute(code), envir = parent.frame())
}
