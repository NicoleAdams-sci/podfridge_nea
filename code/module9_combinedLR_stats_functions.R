# ------------------------------------------------------------------------------
# Module 9: Post-simulation analysis functions on combined LR
# ------------------------------------------------------------------------------
#
# Dependencies:
#   - LR_kinship_utility_functions.R (direct)
#
# Direct function calls:
#   - calculate_summary_stats()
#   - calculate_ratio_stats()
#   - calculate_cutoffs()
#   - calculate_proportions_exceeding_cutoffs()
#
# ------------------------------------------------------------------------------

library(dplyr)
library(tibble)
library(purrr)

# Calculate summary statistics for LRs
calculate_summary_stats <- function(combined_lrs) {
  summary_stats <- combined_lrs |>
    group_by(known_relationship, population, loci_set, tested_population, tested_relationship, is_correct_pop) |>
    summarize(
      n = n(),
      mean_LR = mean(combined_LR, na.rm = TRUE),
      median_LR = median(combined_LR, na.rm = TRUE),
      sd_LR = sd(combined_LR, na.rm = TRUE),
      min_LR = min(combined_LR, na.rm = TRUE),
      max_LR = max(combined_LR, na.rm = TRUE),
      lower_95 = quantile(combined_LR, 0.025, na.rm = TRUE),
      upper_95 = quantile(combined_LR, 0.975, na.rm = TRUE),
      .groups = 'drop'
    ) |>
    ungroup()
  
  return(summary_stats)
}

calculate_ratio_stats <- function(all_combined) {
  all_combined <- as.data.table(all_combined)
  
  # 1. Define Correct Population LRs (LR_C)
  # Ensure the correct LR is calculated using the true population for the hypothesis tested.
  combined_lrs_correct <- all_combined[
    is_correct_pop == TRUE & (tested_population == population),
    .(
      batch_id,
      pair_id,
      population,
      known_relationship,
      tested_relationship,
      loci_set,
      correct_LR = combined_LR
    )
  ]
  
  # 2. Define Wrong Population LRs (LR_W)
  # Ensure the wrong LR is calculated using a population other than the true one.
  combined_lrs_wrong <- all_combined[
    is_correct_pop == FALSE & (tested_population != population),
    .(
      batch_id,
      pair_id,
      population,
      known_relationship,
      tested_relationship,
      loci_set,
      tested_population, # The frequency file used
      wrong_LR = combined_LR
    )
  ]
  
  # 3. De-duplicate the correct set to prevent join inflation (if duplicates exist)
  combined_lrs_correct_unique <- unique(combined_lrs_correct,
                                        by = c("batch_id", "pair_id", "population",
                                               "known_relationship", "tested_relationship", "loci_set"))
  
  # 4. Merge LR_W and LR_C
  # Merge key ensures matching of the unique simulated pair and the specific hypothesis test.
  combined_lrs_ratio <- merge(
    combined_lrs_wrong,
    combined_lrs_correct_unique,
    by = c(
      "batch_id",
      "pair_id",
      "population",
      "known_relationship",
      "loci_set",
      "tested_relationship"
    )
  )
  
  # 5. Calculate the ratio (data.table syntax)
  combined_lrs_ratio[, ratio := wrong_LR / correct_LR]
  
  # 6. Calculate summary statistics (dplyr syntax)
  # Grouping by all experimental factors and the 'wrong' population used.
  ratio_summary <- combined_lrs_ratio %>%
    group_by(
      population,              # True population of the individuals
      known_relationship,      # True relationship of the individuals
      tested_relationship,     # Relationship hypothesis being tested (H1)
      loci_set,                # Loci used for the calculation
      tested_population        # Wrong frequency file used (H2)
    ) %>%
    summarize(
      n = n(),
      mean_ratio = mean(ratio, na.rm = TRUE),
      median_ratio = median(ratio, na.rm = TRUE),
      sd_ratio = sd(ratio, na.rm = TRUE),
      min_ratio = min(ratio, na.rm = TRUE),
      max_ratio = max(ratio, na.rm = TRUE),
      lower_95 = quantile(ratio, 0.025, na.rm = TRUE),
      upper_95 = quantile(ratio, 0.975, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    ungroup()
  
  return(list(ratio_summary = as.data.frame(ratio_summary),
              combined_lrs_ratio = as.data.frame(combined_lrs_ratio)))
}

# Function to calculate cut-off values for 1%, 0.1%, and 0.01% FPR
calculate_cutoffs <- function(all_combined, fp_rates) {
  
  # Note: The input dataframe name is changed to 'all_combined' for consistency
  cutoffs <- all_combined %>%
    # 1. Filter: Identify the distribution of LRs for true "unrelated" pairs
    #    The 'tested_relationship' must be the one corresponding to the LR value (e.g., Parent/Child)
    filter(known_relationship == "unrelated") %>% 
    
    # 2. Group: Calculate separate cutoffs for each tested frequency file and loci set
    group_by(loci_set, tested_population, tested_relationship) %>%
    
    summarize(
      # A fixed, arbitrary cutoff is often included for comparison
      fixed_cutoff = 1.00,
      
      # Calculate cutoffs based on False Positive Rates (FPRs)
      # The quantile function finds the LR value where a given percentage (1 - FPR) of the data falls below it.
      # fp_rates is assumed to be a vector like c(1, 0.1, 0.01)
      
      # Cutoff for 1% FPR (e.g., fp_rates[1]=1)
      cutoff_1 = quantile(combined_LR, probs = 1 - fp_rates[1] / 100, na.rm = TRUE),
      
      # Cutoff for 0.1% FPR (e.g., fp_rates[2]=0.1)
      cutoff_0_1 = quantile(combined_LR, probs = 1 - fp_rates[2] / 100, na.rm = TRUE),
      
      # Cutoff for 0.01% FPR (e.g., fp_rates[3]=0.01)
      cutoff_0_01 = quantile(combined_LR, probs = 1 - fp_rates[3] / 100, na.rm = TRUE),
      
      # Count the number of unrelated pairs used for the calculation
      n_unrelated = n(),
      
      .groups = 'drop'
    )
  
  return(cutoffs)
}

# Function to calculate proportions exceeding cutoffs
calculate_proportions_exceeding_cutoffs <- function(all_combined, cutoffs) {
  
  # 1. Join with cutoffs
  # Key columns are loci_set, tested_population, and tested_relationship (the hypothesis)
  df_with_cutoffs <- left_join(
    # The 'all_combined' data contains the LRs to check
    all_combined,
    # The 'cutoffs' table contains the thresholds
    cutoffs,
    by = c("loci_set", "tested_population", "tested_relationship") # Must join on all grouping variables from cutoffs
  )
  
  # 2. Check which LRs exceed the calculated cutoffs
  df_with_cutoffs <- df_with_cutoffs %>%
    mutate(
      exceeds_fixed_cutoff = combined_LR > fixed_cutoff,
      exceeds_cutoff_1     = combined_LR > cutoff_1,
      exceeds_cutoff_0_1   = combined_LR > cutoff_0_1,
      exceeds_cutoff_0_01  = combined_LR > cutoff_0_01
    )
  
  # 3. Calculate proportions exceeding the cutoffs
  proportions_exceeding <- df_with_cutoffs %>%
    group_by(
      population,              # True population of individuals
      known_relationship,      # True relationship of individuals
      tested_relationship,     # Relationship hypothesis being tested
      loci_set,                # Loci used for calculation
      tested_population,       # Frequency file used for calculation
      is_correct_pop           # Whether the frequency file matches the true population
    ) %>%
    summarize(
      # Calculate proportion by summing TRUE (1) vs. total n()
      proportion_exceeding_fixed = sum(exceeds_fixed_cutoff, na.rm = TRUE) / n(),
      proportion_exceeding_1     = sum(exceeds_cutoff_1, na.rm = TRUE) / n(),
      proportion_exceeding_0_1   = sum(exceeds_cutoff_0_1, na.rm = TRUE) / n(),
      proportion_exceeding_0_01  = sum(exceeds_cutoff_0_01, na.rm = TRUE) / n(),
      n_related = n(),
      .groups = 'drop'
    ) %>%
    # Filter out 'unrelated' pairs as they were used to define the cutoffs (FPR)
    # This result table focuses on the True Positive Rate (TPR) for related pairs.
    filter(known_relationship != "unrelated")
  
  return(proportions_exceeding)
}