##### Test Module flow #####
library(dplyr)
library(data.table)

# Source the utility functions to get FALLBACK_FREQ and other shared constants
source("code/LR_kinship_utility_functions.R")
source("code/module1_allele_simulator.R")
source("code/module2_STR_profile_simulator.R")

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
