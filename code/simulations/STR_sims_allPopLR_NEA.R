# Load Required Libraries
suppressMessages(suppressWarnings({
  library(tidyverse)
  library(furrr)
  library(progressr)
  library(data.table)
  library(ggplot2)
  library(future)
  library(parallel)
}))

# Usage: gets submitted by script_NEA.sh (Rscript code/STR_sims_allPopLR_NEA.R [RELATED] [UNRELATED])

# Set up cluster
ncores <- as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK', unset = 4))

cl <- makeCluster(ncores)
plan(cluster, workers = cl)

# Ensure the cluster is stopped when the script exits
on.exit(parallel::stopCluster(cl))

# Helper function for logging
log_message <- function(message) {
  cat(paste0("[", Sys.time(), "] ", message, "\n"))
}

# Helper function to log function timings
timing_log <- list()

log_function_time <- function(func, name, ...) {
  start_time <- Sys.time()
  result <- func(...)
  end_time <- Sys.time()
  duration <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  if (!name %in% names(timing_log)) {
    timing_log[[name]] <- list(total = 0, count = 0, min = Inf, max = -Inf, times = c())
  }
  
  timing_log[[name]]$total <- timing_log[[name]]$total + duration
  timing_log[[name]]$count <- timing_log[[name]]$count + 1
  timing_log[[name]]$min <- min(timing_log[[name]]$min, duration)
  timing_log[[name]]$max <- max(timing_log[[name]]$max, duration)
  timing_log[[name]]$times <- c(timing_log[[name]]$times, duration)
  
  return(result)
}

# Read Command-Line Arguments
args <- commandArgs(trailingOnly = TRUE)
n_sims_related <- as.numeric(args[1])
n_sims_unrelated <- as.numeric(args[2])

# Get SLURM array job ID and task ID ** changed here ***
array_job_id <- Sys.getenv("SLURM_ARRAY_JOB_ID", unset = "single")
task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
job_id <- Sys.getenv("SLURM_JOB_ID")

# Create output folder with timestamp
timestamp <- format(Sys.time(), "%Y%m%d")
# ** changed here ***
output_dir <- file.path("output", paste0("simulation_", timestamp, "_array", task_id))
dir.create(output_dir, recursive = TRUE)

# ** changed here ***
output_file <- file.path(output_dir, paste0("sim_processed_genotypes_task", task_id, ".csv"))
summary_output_file <- file.path(output_dir, paste0("sim_summary_genotypes_task", task_id, ".csv"))
timing_log_file <- file.path(output_dir, paste0("timing_log_task", task_id, ".csv"))

# Log the start of the process
log_message(paste("Starting simulation setup and processing for array task", task_id))

# Load Allele Frequencies Data
log_message("Loading allele frequencies data...")
allele_freq_time <- system.time({
  df_allelefreq <- fread("data/df_allelefreq_combined.csv")
  df_allelefreq <- df_allelefreq[population != "all"] # Filter out "all" population
  df_allelefreq$frequency <- ifelse(df_allelefreq$frequency==0, 5/(2*1036), df_allelefreq$frequency) # add minimum allele freq based on forensic best practices
  df_allelefreq[, allele := as.character(allele)]
})
log_message(paste("Loaded allele frequencies data in", allele_freq_time["elapsed"], "seconds."))

# Extract unique loci
log_message("Extracting unique loci...")
loci_list <- df_allelefreq |>
  pull(marker) |>
  unique()

# Load Core Loci Data
log_message("Loading core loci data...")
core_loci_time <- system.time({
  core_loci <- fread("data/core_CODIS_loci.csv")
  columns <- c("core_13", "identifiler_15", "expanded_20", "supplementary")
  loci_lists <- lapply(columns, function(col) {
    core_loci |>
      filter(get(col) == 1) |>
      pull(locus)
  })
  names(loci_lists) <- columns
  loci_lists$autosomal_29 <- loci_list
})
log_message(paste("Loaded core loci data in", core_loci_time["elapsed"], "seconds."))

# Define Kinship Matrix
kinship_matrix <- tibble(
  relationship_type = factor(
    c("parent_child", "full_siblings", "half_siblings", "cousins", "second_cousins", "unrelated"),
    levels = c("parent_child", "full_siblings", "half_siblings", "cousins", "second_cousins", "unrelated")
  ),
  k0 = c(0, 1/4, 1/2, 7/8, 15/16, 1),
  k1 = c(1, 1/2, 1/2, 1/8, 1/16, 0),
  k2 = c(0, 1/4, 0, 0, 0, 0)
)

# Define Populations
population_labels <- tibble(
  population = factor(
    c("AfAm", "Cauc", "Hispanic", "Asian"),
    levels = c("AfAm", "Cauc", "Hispanic", "Asian")
  ),
  label = c("African American", "Caucasian", "Hispanic", "Asian")
)
populations_list <- levels(population_labels$population)

# Functions
calculate_likelihood_ratio <- function(shared_alleles, genotype_match = NULL, pA = NULL, pB = NULL, k0, k1, k2) {
  if (shared_alleles == 0) {
    LR <- k0
    return(LR)
  }
  if (shared_alleles == 1) {
    if (genotype_match == "AA-AA") {
      Rxp <- pA
    } else if (genotype_match == "AA-AB" | genotype_match == "AB-AA") {
      Rxp <- 2 * pA
    } else if (genotype_match == "AB-AC") {
      Rxp <- 4 * pA
    } else if (genotype_match == "AB-AB") {
      Rxp <- (4 * (pA * pB)) / (pA + pB)
    } else {
      stop("Invalid genotype match for 1 shared allele.")
    }
    LR <- k0 + (k1 / Rxp)
    return(LR)
  }
  if (shared_alleles == 2) {
    if (genotype_match == "AA-AA") {
      Rxp <- pA
      Rxu <- pA^2
    } else if (genotype_match == "AB-AB") {
      Rxp <- (4 * pA * pB) / (pA + pB)
      Rxu <- 2 * pA * pB
    } else {
      stop("Invalid genotype match for 2 shared alleles.")
    }
    LR <- k0 + (k1 / Rxp) + (k2 / Rxu)
    return(LR)
  }
}

generate_simulation_setup <- function(kinship_matrix, population_list, num_related, num_unrelated) {
  simulation_setup <- data.frame(
    population = character(),
    relationship_type = character(),
    num_simulations = integer(),
    stringsAsFactors = FALSE
  )
  for (population in population_list) {
    for (relationship in kinship_matrix$relationship_type) {
      num_simulations <- ifelse(relationship == "unrelated", num_unrelated, num_related)
      simulation_setup <- rbind(simulation_setup, data.frame(
        population = population,
        relationship_type = relationship,
        num_simulations = num_simulations
      ))
    }
  }
  return(simulation_setup)
}

# Initialize individuals pair with LR columns for ALL populations - to test known v tested
initialize_individuals_pair <- function(population, relationship_type, sim_id, loci_list, population_list) {
  num_loci <- length(loci_list)
  individuals_genotypes <- data.table(
    population = rep(population, num_loci),
    relationship_type = rep(relationship_type, num_loci),
    sim_id = rep(sim_id, num_loci),
    locus = loci_list,
    ind1_allele1 = character(num_loci),
    ind1_allele2 = character(num_loci),
    ind2_allele1 = character(num_loci),
    ind2_allele2 = character(num_loci),
    shared_alleles = integer(num_loci),
    genotype_match = character(num_loci)
  )
  
  # Add columns for LRs calculated using ALL populations' allele frequencies
  for (pop in population_list) {
    individuals_genotypes[, paste0("LR_", pop) := numeric(num_loci)]
  }
  
  return(individuals_genotypes)
}

simulate_genotypes <- function(row, df_allelefreq, kinship_matrix) {
  population <- row$population
  locus <- row$locus
  relationship <- row$relationship_type
  
  allele_freqs <- df_allelefreq |>
    filter(population == !!population, marker == !!locus, frequency > 0)
  
  if (nrow(allele_freqs) == 0) {
    stop(paste("No valid alleles found for population", population, "and locus", locus))
  }
  
  alleles <- allele_freqs$allele
  frequencies <- allele_freqs$frequency
  
  frequencies <- round(frequencies / sum(frequencies), 6)
  valid_indices <- frequencies > 0
  alleles <- alleles[valid_indices]
  frequencies <- frequencies[valid_indices]
  
  ind1_alleles <- sample(alleles, size = 2, replace = TRUE, prob = frequencies)
  
  kinship_coeffs <- kinship_matrix[kinship_matrix$relationship_type == relationship, ]
  relationship_choice <- sample(c('none', 'one', 'both'), size = 1, prob = c(kinship_coeffs$k0, kinship_coeffs$k1, kinship_coeffs$k2))
  
  if (relationship_choice == 'none') {
    ind2_alleles <- sample(alleles, size = 2, replace = TRUE, prob = frequencies)
  } else if (relationship_choice == 'one') {
    shared_allele <- sample(ind1_alleles, size = 1)
    non_shared_allele <- sample(alleles, size = 1, replace = TRUE, prob = frequencies)
    if (runif(1) > 0.5) {
      ind2_alleles <- c(shared_allele, non_shared_allele)
    } else {
      ind2_alleles <- c(non_shared_allele, shared_allele)
    }
  } else if (relationship_choice == 'both') {
    ind2_alleles <- ind1_alleles
  }
  
  row$ind1_allele1 <- ind1_alleles[1]
  row$ind1_allele2 <- ind1_alleles[2]
  row$ind2_allele1 <- ind2_alleles[1]
  row$ind2_allele2 <- ind2_alleles[2]
  
  return(row)
}

# Standardized kinship_calculation function that calculates LRs for ALL populations
kinship_calculation <- function(row, allele_frequency_data, kinship_matrix, population_list) {
  alleles_ind1 <- as.character(c(row$ind1_allele1, row$ind1_allele2))
  alleles_ind2 <- as.character(c(row$ind2_allele1, row$ind2_allele2))
  
  shared_alleles_vector <- intersect(alleles_ind1, alleles_ind2)
  unique_alleles_ind1 <- setdiff(alleles_ind1, shared_alleles_vector)
  unique_alleles_ind2 <- setdiff(alleles_ind2, shared_alleles_vector)
  
  allele_map <- list()
  next_label <- 1
  
  for (allele in shared_alleles_vector) {
    allele_map[[LETTERS[next_label]]] <- allele
    next_label <- next_label + 1
  }
  
  for (allele in unique_alleles_ind1) {
    if (!(allele %in% allele_map)) {
      allele_map[[LETTERS[next_label]]] <- allele
      next_label <- next_label + 1
    }
  }
  
  for (allele in unique_alleles_ind2) {
    if (!(allele %in% allele_map)) {
      allele_map[[LETTERS[next_label]]] <- allele
      next_label <- next_label + 1
    }
  }
  
  allele_map <- unlist(allele_map)
  labeled_alleles_ind1 <- sapply(as.character(alleles_ind1), function(x) names(allele_map)[which(allele_map == x)])
  labeled_alleles_ind2 <- sapply(as.character(alleles_ind2), function(x) names(allele_map)[which(allele_map == x)])
  
  shared_alleles <- length(shared_alleles_vector)
  genotype_ind1 <- paste(sort(labeled_alleles_ind1), collapse = "")
  genotype_ind2 <- paste(sort(labeled_alleles_ind2), collapse = "")
  genotype_match <- paste(genotype_ind1, genotype_ind2, sep = "-")
  
  row$shared_alleles <- shared_alleles
  row$genotype_match <- genotype_match
  
  # Calculate LR for ALL populations (no special case for correct population)
  A_allele <- ifelse("A" %in% names(allele_map), allele_map[["A"]], NA)
  B_allele <- ifelse("B" %in% names(allele_map), allele_map[["B"]], NA)
  
  # Get kinship coefficients
  k_values <- kinship_matrix[kinship_matrix$relationship_type == row$relationship_type, ]
  
  # Calculate LR for each population
  for (pop in population_list) {
    pop_freqs <- dplyr::filter(allele_frequency_data, population == pop, marker == row$locus)
    
    if (nrow(pop_freqs) > 0) {
      pA <- ifelse(any(pop_freqs$allele == A_allele), pop_freqs$frequency[pop_freqs$allele == A_allele], NA)
      pB <- ifelse(any(pop_freqs$allele == B_allele), pop_freqs$frequency[pop_freqs$allele == B_allele], NA)
      
      if (!is.na(pA)) {
        if (is.na(pB) && length(shared_alleles_vector) > 1) {
          # If we need pB but it's missing, set LR to NA
          row[[paste0("LR_", pop)]] <- NA
        } else {
          row[[paste0("LR_", pop)]] <- calculate_likelihood_ratio(shared_alleles, genotype_match, pA, pB, k_values$k0, k_values$k1, k_values$k2)
        }
      } else {
        row[[paste0("LR_", pop)]] <- NA
      }
    } else {
      row[[paste0("LR_", pop)]] <- NA
    }
  }
  
  return(row)
}


process_loci <- function(row, allele_frequency_data, kinship_matrix, population_list) {
  simulated_row <- log_function_time(simulate_genotypes, "simulate_genotypes", row, allele_frequency_data, kinship_matrix)
  final_row <- log_function_time(kinship_calculation, "kinship_calculation", simulated_row, allele_frequency_data, kinship_matrix, population_list)
  return(final_row)
}


process_individuals_genotypes <- function(individuals_genotypes, df_allelefreq, kinship_matrix, population_list) {
  final_individuals_genotypes <- individuals_genotypes |>
    future_pmap(~ log_function_time(process_loci, "process_loci", list(...), df_allelefreq, kinship_matrix, population_list), seed = TRUE) |>
    bind_rows()
  return(final_individuals_genotypes)
}

# Standardized approach for combined LRs
# calculate_combined_lrs <- function(final_results, loci_lists, population_list) {
#   final_results <- as.data.table(final_results)
#   
#   # Create list of expressions for calculating products
#   calc_expressions <- list()
#   
#   # For each population LR
#   for (pop in population_list) {
#     col_name <- paste0("LR_", pop)
#     calc_expressions[[paste0("core_13_", pop)]] <- parse(text = paste0('prod(', col_name, '[locus %in% loci_lists$core_13], na.rm = TRUE)'))
#     calc_expressions[[paste0("identifiler_15_", pop)]] <- parse(text = paste0('prod(', col_name, '[locus %in% loci_lists$identifiler_15], na.rm = TRUE)'))
#     calc_expressions[[paste0("expanded_20_", pop)]] <- parse(text = paste0('prod(', col_name, '[locus %in% loci_lists$expanded_20], na.rm = TRUE)'))
#     calc_expressions[[paste0("supplementary_", pop)]] <- parse(text = paste0('prod(', col_name, '[locus %in% loci_lists$supplementary], na.rm = TRUE)'))
#     calc_expressions[[paste0("autosomal_29_", pop)]] <- parse(text = paste0('prod(', col_name, '[locus %in% loci_lists$autosomal_29], na.rm = TRUE)'))
#   }
#   
#   # Execute all calculations
#   combined_lrs <- final_results[, c(as.list(eval(calc_expressions))), by = .(population, relationship_type, sim_id)]
#   
#   # Melt the data
#   measure_vars <- names(combined_lrs)[!(names(combined_lrs) %in% c("population", "relationship_type", "sim_id"))]
#   combined_lrs <- melt(combined_lrs,
#                        id.vars = c("population", "relationship_type", "sim_id"),
#                        measure.vars = measure_vars,
#                        variable.name = "loci_set", 
#                        value.name = "LR")
#   
#   # Parse the loci_set to extract base loci set and frequency population
#   combined_lrs[, c("base_loci_set", "freq_population") := tstrsplit(as.character(loci_set), "_", fixed = TRUE, keep = c(1, 2))]
#   
#   # Add a column to identify if this is the correct population match
#   combined_lrs[, is_correct_pop := (population == freq_population)]
#   
#   return(combined_lrs)
# }

# Updated plotting function for standardized approach
# plot_and_save_results <- function(combined_lrs) {
#   log_message("Starting to plot LR distributions...")
#   
#   # Ensure factor levels are set correctly for plotting
#   combined_lrs$relationship_type <- factor(combined_lrs$relationship_type, levels = c("parent_child", "full_siblings", "half_siblings", "cousins", "second_cousins", "unrelated"))
#   
#   # Plot by assigned population vs frequency population
#   ggplot(combined_lrs, aes(x = relationship_type, y = LR, fill = freq_population, color = freq_population)) +
#     geom_boxplot(position = position_dodge(width = 0.9)) +
#     facet_grid(population ~ base_loci_set, scales = "free_y") +
#     labs(
#       title = "LR Distributions: Assigned Population vs. Frequency Population",
#       x = "Relationship Type",
#       y = "LR (log scale)",
#       fill = "Frequency Population",
#       color = "Frequency Population"
#     ) +
#     theme_minimal() +
#     theme(
#       axis.text.x = element_text(angle = 45, hjust = 1)
#     ) +
#     scale_y_log10()
#   
#   ggsave(file.path(output_dir, "sim_log_lr_panel_plot_by_freq_pop.png"), width = 16, height = 12)
#   
#   # Special highlight for correct population
#   ggplot(combined_lrs, aes(x = relationship_type, y = LR, fill = is_correct_pop)) +
#     geom_boxplot() +
#     facet_grid(population ~ base_loci_set, scales = "free_y") +
#     labs(
#       title = "LR Distributions: Correct vs. Incorrect Population Frequencies",
#       x = "Relationship Type",
#       y = "LR (log scale)",
#       fill = "Correct Population"
#     ) +
#     theme_minimal() +
#     theme(
#       axis.text.x = element_text(angle = 45, hjust = 1)
#     ) +
#     scale_y_log10() +
#     scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "forestgreen"))
#   
#   ggsave(file.path(output_dir, "sim_log_lr_correct_vs_incorrect.png"), width = 16, height = 12)
#   
#   log_message("LR distributions plot saved.")
#   
#   # Calculate and plot summary statistics
#   summary_stats <- combined_lrs |>
#     group_by(relationship_type, population, base_loci_set, freq_population, is_correct_pop) |>
#     summarize(
#       mean_LR = mean(LR, na.rm = TRUE),
#       median_LR = median(LR, na.rm = TRUE),
#       lower_95 = quantile(LR, 0.025, na.rm = TRUE),
#       upper_95 = quantile(LR, 0.975, na.rm = TRUE),
#       .groups = 'drop'
#     ) |>
#     ungroup()
#   
#   # Save summary statistics
#   fwrite(summary_stats, file.path(output_dir, "sim_summary_stats_by_freq_pop.csv"))
#   
#   # Plot ratio of wrong vs correct LR
#   # Create a wide format with correct and wrong LRs side by side
#   combined_lrs_wide <- dcast(
#     combined_lrs[, .(population, relationship_type, sim_id, base_loci_set, freq_population, LR, is_correct_pop)],
#     population + relationship_type + sim_id + base_loci_set ~ is_correct_pop,
#     value.var = "LR"
#   )
#   
#   # Rename columns for clarity
#   setnames(combined_lrs_wide, c("TRUE", "FALSE"), c("correct_LR", "wrong_LR"))
#   
#   # Reshape to get one row per wrong population
#   combined_lrs_wrong <- combined_lrs[is_correct_pop == FALSE]
#   combined_lrs_correct <- combined_lrs[is_correct_pop == TRUE, .(population, relationship_type, sim_id, base_loci_set, correct_LR = LR)]
#   
#   combined_lrs_ratio <- merge(
#     combined_lrs_wrong,
#     combined_lrs_correct,
#     by = c("population", "relationship_type", "sim_id", "base_loci_set")
#   )
#   
#   # Calculate ratio of wrong to correct LR
#   combined_lrs_ratio[, ratio := LR / correct_LR]
#   
#   # Plot ratio distribution
#   ggplot(combined_lrs_ratio, aes(x = relationship_type, y = ratio, fill = freq_population)) +
#     geom_boxplot() +
#     facet_grid(population ~ base_loci_set, scales = "free_y") +
#     labs(
#       title = "Ratio of Wrong Population LR to Correct Population LR",
#       x = "Relationship Type",
#       y = "Wrong LR / Correct LR (log scale)",
#       fill = "Wrong Population"
#     ) +
#     theme_minimal() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     scale_y_log10()
#   
#   ggsave(file.path(output_dir, "sim_lr_ratio_plot.png"), width = 16, height = 12)
#   
#   # Save ratio summary statistics
#   ratio_summary <- combined_lrs_ratio |>
#     group_by(relationship_type, population, base_loci_set, freq_population) |>
#     summarize(
#       mean_ratio = mean(ratio, na.rm = TRUE),
#       median_ratio = median(ratio, na.rm = TRUE),
#       lower_95 = quantile(ratio, 0.025, na.rm = TRUE),
#       upper_95 = quantile(ratio, 0.975, na.rm = TRUE),
#       .groups = 'drop'
#     ) |>
#     ungroup()
#   
#   fwrite(ratio_summary, file.path(output_dir, "sim_lr_ratio_summary.csv"))
# }

# Updated cutoffs calculation for standardized approach
# calculate_cutoffs <- function(input_df, fp_rates) {
#   cutoffs <- input_df |>
#     filter(relationship_type == "unrelated") |>
#     group_by(base_loci_set, freq_population) |>
#     summarize(
#       fixed_cutoff = 1.00,
#       cutoff_1 = quantile(LR, probs = 1 - fp_rates[1] / 100, na.rm = TRUE),
#       cutoff_0_1 = quantile(LR, probs = 1 - fp_rates[2] / 100, na.rm = TRUE),
#       cutoff_0_01 = quantile(LR, probs = 1 - fp_rates[3] / 100, na.rm = TRUE),
#       n_unrelated = n(),
#       .groups = 'drop'
#     )
#   return(cutoffs)
# }

# Updated proportions exceeding cutoffs calculation for standardized approach
# calculate_proportions_exceeding_cutoffs <- function(input_df, cutoffs) {
#   # Join with cutoffs  
#   df_with_cutoffs <- left_join(
#     input_df, 
#     cutoffs, 
#     by = c("base_loci_set", "freq_population")
#   )
#   
#   df_with_cutoffs <- df_with_cutoffs |>
#     mutate(
#       exceeds_fixed_cutoff = LR > fixed_cutoff,
#       exceeds_cutoff_1 = LR > cutoff_1,
#       exceeds_cutoff_0_1 = LR > cutoff_0_1,
#       exceeds_cutoff_0_01 = LR > cutoff_0_01
#     )
#   
#   proportions_exceeding <- df_with_cutoffs |>
#     group_by(population, relationship_type, base_loci_set, freq_population, is_correct_pop) |>
#     summarize(
#       proportion_exceeding_fixed = sum(exceeds_fixed_cutoff, na.rm = TRUE) / n(),
#       proportion_exceeding_1 = sum(exceeds_cutoff_1, na.rm = TRUE) / n(),
#       proportion_exceeding_0_1 = sum(exceeds_cutoff_0_1, na.rm = TRUE) / n(),
#       proportion_exceeding_0_01 = sum(exceeds_cutoff_0_01, na.rm = TRUE) / n(),
#       n_related = n(),
#       .groups = 'drop'
#     ) |>
#     filter(relationship_type != "unrelated")
#   
#   return(proportions_exceeding)
# }

# Updated plot proportions exceeding cutoffs for standardized approach
# plot_proportions_exceeding_cutoffs <- function(proportions_exceeding_cutoffs) {
#   log_message("Starting to plot proportions exceeding cutoffs...")
#   
#   # Ensure factor levels are set correctly for plotting
#   proportions_exceeding_cutoffs$relationship_type <- factor(proportions_exceeding_cutoffs$relationship_type, levels = c("parent_child", "full_siblings", "half_siblings", "cousins", "second_cousins", "unrelated"))
#   proportions_exceeding_cutoffs$population <- factor(proportions_exceeding_cutoffs$population, levels = c("AfAm", "Cauc", "Hispanic", "Asian"))
#   
#   # Create long format data for plotting
#   proportions_long <- proportions_exceeding_cutoffs |>
#     pivot_longer(cols = starts_with("proportion_exceeding"),
#                  names_to = "Cutoff_Type", values_to = "Proportion",
#                  names_prefix = "proportion_exceeding_")
#   
#   proportions_long$Cutoff_Type <- factor(proportions_long$Cutoff_Type, levels = c("fixed", "1", "0_1", "0_01"),
#                                          labels = c("Fixed Cutoff (1.00)", "1% FPR", "0.1% FPR", "0.01% FPR"))
#   
#   # Plot proportions with frequency population
#   ggplot(proportions_long, aes(x = relationship_type, y = Proportion, fill = freq_population)) +
#     geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
#     facet_grid(population ~ Cutoff_Type + base_loci_set, scales = "fixed") +
#     labs(
#       title = "Proportions Exceeding Likelihood Cut-offs: By Population and Frequency Source",
#       x = "Relationship Type",
#       y = "Proportion Exceeding Cut-off",
#       fill = "Frequency Population"
#     ) +
#     theme_minimal() +
#     theme(
#       axis.text.x = element_text(angle = 45, hjust = 1)
#     )
#   
#   ggsave(file.path(output_dir, "sim_proportions_exceeding_cutoffs_by_freq_pop.png"), width = 16, height = 12)
#   
#   # Plot proportions with correct vs incorrect frequency
#   ggplot(proportions_long, aes(x = relationship_type, y = Proportion, fill = is_correct_pop)) +
#     geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
#     facet_grid(population ~ Cutoff_Type + base_loci_set, scales = "fixed") +
#     labs(
#       title = "Proportions Exceeding Likelihood Cut-offs: Correct vs. Incorrect Population",
#       x = "Relationship Type",
#       y = "Proportion Exceeding Cut-off",
#       fill = "Correct Population"
#     ) +
#     theme_minimal() +
#     theme(
#       axis.text.x = element_text(angle = 45, hjust = 1)
#     ) +
#     scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "forestgreen"))
#   
#   ggsave(file.path(output_dir, "sim_proportions_exceeding_cutoffs_correct_vs_incorrect.png"), width = 16, height = 12)
#   
#   log_message("Proportions exceeding cutoffs plots saved.")
# }

# Main process function
process_simulation_setup <- function(simulation_setup, df_allelefreq, kinship_matrix, loci_list, loci_lists, output_file, summary_output_file) {
  log_message("Processing simulation setup...")
  process_time <- system.time({
    task_id_num <- as.numeric(task_id) #Get the task_id to use as an offset for sim_id for array
    final_results <- simulation_setup |>
      future_pmap_dfr(function(population, relationship_type, num_simulations) {
        purrr::map_dfr(1:num_simulations, function(sim_id) {
          global_sim_id <- sim_id + (task_id_num - 1) * num_simulations # Modify sim_id to be unique across arrays
          individuals_genotypes <- initialize_individuals_pair(population, relationship_type, global_sim_id, loci_list, populations_list)
          processed_genotypes <- log_function_time(process_individuals_genotypes, "process_individuals_genotypes", individuals_genotypes, df_allelefreq, kinship_matrix, populations_list)
          return(processed_genotypes)
        })
      })
    
    if ("seed" %in% colnames(final_results)) {
      final_results <- final_results |> select(-seed)
    }
    
    # Save the raw genotype results with all population LRs
    fwrite(final_results, output_file)
    
    # Calculate combined LRs and save
    #combined_lrs <- log_function_time(calculate_combined_lrs, "calculate_combined_lrs", final_results, loci_lists, populations_list)
    #fwrite(combined_lrs, summary_output_file)
    
    # Create plots and save
    #log_function_time(plot_and_save_results, "plot_and_save_results", combined_lrs)
    
    # Calculate and save cutoffs - now using standardized approach
    #cutoffs <- log_function_time(calculate_cutoffs, "calculate_cutoffs", combined_lrs, c(1, 0.1, 0.01))
    #fwrite(cutoffs, file.path(output_dir, "sim_cutoffs.csv"))
    
    # Calculate proportions exceeding cutoffs - now using standardized approach
    #proportions_exceeding_cutoffs <- log_function_time(calculate_proportions_exceeding_cutoffs, "calculate_proportions_exceeding_cutoffs", combined_lrs, cutoffs)
    #fwrite(proportions_exceeding_cutoffs, file.path(output_dir, "sim_proportions_exceeding_cutoffs.csv"))
    
    # Plot proportions exceeding cutoffs - now using standardized approach
    #log_function_time(plot_proportions_exceeding_cutoffs, "plot_proportions_exceeding_cutoffs", proportions_exceeding_cutoffs)
  })
  log_message(paste("Processing completed in", process_time["elapsed"], "seconds."))
}

# Execute Simulation Setup and Processing
simulation_setup <- log_function_time(generate_simulation_setup, "generate_simulation_setup", kinship_matrix, populations_list, n_sims_related, n_sims_unrelated)
log_function_time(process_simulation_setup, "process_simulation_setup", simulation_setup, df_allelefreq, kinship_matrix, loci_list, loci_lists, output_file, summary_output_file)
log_message("Simulation setup and processing completed.")

# Save timing log to CSV
timing_log_df <- tibble(
  function_name = names(timing_log),
  total_time = sapply(timing_log, function(x) x$total),
  count = sapply(timing_log, function(x) x$count),
  min_time = sapply(timing_log, function(x) x$min),
  max_time = sapply(timing_log, function(x) x$max),
  avg_time = sapply(timing_log, function(x) x$total / x$count)
)

write_csv(timing_log_df, timing_log_file)