##### Test Creating Pairs #####
# USAGE: module load Rtidyverse; Rscript code/sim_pairs_test.R <POP> <RELATIONSHIP> <No. pairs>
# Example usage: module load Rtidyverse; Rscript code/sim_pairs_test.R AfAm full_siblings 3


# Load in arguments
args <- commandArgs(trailingOnly = TRUE)
POP <- args[1]
SHIP <- args[2]
N_PAIRS <- ifelse(length(args) >= 3, as.numeric(args[3]), 1)  # Default N=1 if not provided

# Load packages
suppressPackageStartupMessages({
	library(dplyr)
	library(data.table)
	library(tictoc)
})


# Source the utility functions to get FALLBACK_FREQ and other shared constants
source("code/LR_kinship_utility_functions.R")
source("code/module1_allele_simulator.R")
source("code/module2_STR_profile_simulator.R")
source("code/module3_related_individual_simulator.R")
source("code/module4_single_locus_LR.R")
# source("code/module5_combined_LR.R")
# source("code/module6_single_combo_pair_generator.R")

#### Load in necessary datasets ####
# Load kinship matrix
kinship_matrix <- fread("data/kinship_coefficients.csv")


#### Simulate Pairs ####
# Use Module 3 - Simulate multiple pairs
print(paste("Simulating", N_PAIRS, "pairs for", POP, SHIP))
tic("simulate_individual_pair")

if (N_PAIRS == 1) {
  # Use single pair function for consistency with original
  pairs <- simulate_individual_pair(loci_list, POP, SHIP, df_allelefreq)
} else {
  # Use multiple pairs function
  pairs <- simulate_multiple_pairs(N_PAIRS, loci_list, POP, SHIP, df_allelefreq)
}

toc()

filenam <- paste("pairs", POP, SHIP, paste0("n", N_PAIRS), format(Sys.time(), "%Y%m%d"), sep="_")

write.csv(pairs, file=paste0("output/",filenam, ".csv"), quote=FALSE, row.names = FALSE)
