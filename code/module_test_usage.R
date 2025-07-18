##### Test Module flow #####
library(dplyr)
library(data.table)

setwd("~/Documents/lasisiLab/podfridge_nea/")

# Source the utility functions to get FALLBACK_FREQ and other shared constants
source("code/LR_kinship_utility_functions.R")
source("code/module1_allele_simulator.R")
source("code/module2_STR_profile_simulator.R")
source("code/module4_related_individual_simulator.R")
source("code/module5_single_locus_LR.R")
source("code/module6_combined_LR.R")
source("code/module7_single_combo_pair_generator.R")
source("code/module8_single_pop_focal_generator.R")
source("code/module9_unrelated_pool_generator.R")

# Load Allele Frequencies Data
df_allelefreq <- fread("data/df_allelefreq_combined.csv")
#df_allelefreq <- df_allelefreq[population != "all"] # Filter out "all" population
df_allelefreq$frequency <- ifelse(df_allelefreq$frequency == 0, FALLBACK_FREQ, df_allelefreq$frequency) # Use shared constant
df_allelefreq[, allele := as.character(allele)]

# Extract unique loci
loci_list <- df_allelefreq |> pull(marker) |> unique()

# Load Core Loci Data
core_loci <- fread("data/core_CODIS_loci.csv")
columns <- c("core_13", "identifiler_15", "expanded_20", "supplementary")
loci_lists <- lapply(columns, function(col) {
  core_loci |>
    filter(get(col) == 1) |>
    pull(locus)
})
names(loci_lists) <- columns
loci_lists$autosomal_29 <- loci_list

# Load kinship matrix
kinship_matrix <- fread("data/kinship_coefficients.csv")

# Use Module 1 - simulate single allele
single_allele <- simulate_allele("D3S1358", "AfAm", df_allelefreq)

# Use Module 2 - simulate STR profile (This could be core_13, identifiler_15, or any custom list)
profile <- simulate_str_profile(loci_lists$core_13, "AfAm", df_allelefreq)

# Use Module 4 -# Simulate multiple pairs
pairs <- simulate_individual_pair(loci_list, "AfAm", "full_siblings", df_allelefreq)

# Use Module 5 - calculate single locus LR 
tested_relationships <- c("parent_child", "full_siblings", "unrelated")
tested_populations <- c("AfAm", "Cauc", "Hispanic", "Asian")

# Calculate single-locus LRs
single_locus_lr <- calculate_single_locus_lr(
  pair_data = pairs,
  tested_relationship = "full_siblings",
  tested_populations = tested_populations,
  allele_frequency_data = df_allelefreq,
  kinship_coefficients = kinship_matrix
)

# Use Module 6 - calculate loci set LR
combined_results <- calculate_combined_lr(single_locus_lr, loci_lists)

# Use Module 7 
# Single combination
result <- generate_single_combo_pairs(
  population = "AfAm",
  relationship = "full_siblings",
  n_pairs = 100,
  loci_list = loci_list,
  allele_frequency_data = df_allelefreq,
  kinship_coefficients = kinship_matrix
)

# Use Module 8
# Single population focal generation
focalresults <- generate_single_pop_focal(
  population = "AfAm",
  n_focal = 100,
  family_structure = "nuclear",
  loci_list = loci_list,
  allele_frequency_data = df_allelefreq,
  kinship_coefficients = kinship_matrix
)

# # Multiple populations
# focalresults2 <- generate_multiple_pop_focal(
#   populations = c("AfAm", "Cauc", "Hispanic", "Asian"),
#   n_focal_per_pop = 50,
#   family_structure = "siblings_only",
#   loci_list = loci_list,
#   allele_frequency_data = df_allelefreq,
#   kinship_coefficients = kinship_matrix
# )

# Use Module 9
# Single population unrelated pool
unrelatedresult <- generate_unrelated_pool(
  population = "AfAm",
  n_unrelated = 1000,
  loci_list = loci_list,
  allele_frequency_data = df_allelefreq
)

# # Multiple populations
# unrelatedresult <- generate_multiple_pop_unrelated(
#   populations = c("AfAm", "Cauc", "Hispanic", "Asian"),
#   n_unrelated_per_pop = 1000,
#   loci_list = loci_list,
#   allele_frequency_data = df_allelefreq
# )
