# PODFRIDGE Simulation Pipeline

## Overview

The PODFRIDGE (POPulation Data For Relationship Inference using DNA Genotype Evidence) pipeline simulates STR (Short Tandem Repeat) profiles for genetically related individual pairs and calculates likelihood ratios (LRs) for kinship analysis. These simulations are essential for understanding how different population assumptions and relationship hypotheses affect forensic DNA kinship testing accuracy.

### Why These Simulations?

1. **Evaluate Population Effects**: Assess how using incorrect population allele frequencies impacts LR calculations
2. **Test Relationship Misspecification**: Quantify errors when testing wrong relationship hypotheses
3. **Optimize Loci Sets**: Compare performance of different STR marker panels (Core 13, Identifiler 15, Expanded 20, etc.)
4. **Establish Decision Thresholds**: Determine appropriate LR cutoffs for different false positive rates

## Directory Structure

```
podfridge_simulations/
├── code/                    # All R scripts and SLURM job scripts
├── data/                    # Reference data files
│   ├── df_allelefreq_combined.csv    # Allele frequencies for 5 populations
│   ├── kinship_coefficients.csv      # IBD sharing probabilities (k0, k1, k2)
│   └── core_CODIS_loci.csv          # Loci set definitions
├── logs/                    # SLURM job logs and error files
└── output/                  # Generated data and results
    ├── pairs_*.csv          # Simulated genotype pairs
    ├── LR/                  # Single-locus likelihood ratios
    ├── combined_LR/         # Multi-locus combined LRs
    └── lr_analysis_*/       # Analysis results and plots
        ├── plots_matched/   # Plots for correct population/relationship
        ├── plots_mismatched/# Plots for incorrect assumptions
        └── plots_exceeding_cutoffs/  # Threshold analysis plots
```

## Module Architecture & Dependencies

### Core Simulation Modules (R Scripts)

```
Foundation:
└── LR_kinship_utility_functions.R
    ├── Provides: FALLBACK_FREQ constant, loci_lists, df_allelefreq
    ├── Helper functions: count_shared_alleles(), label_and_genotype()
    └── Main function: kinship_calculation()

Simulations:
└── module1_allele_simulator.R
    └── module2_STR_profile_simulator.R
        └── module3_related_individual_simulator.R

Analyses:
└── module4_single_locus_LR.R
    └── module5_combined_LR.R

Post-Processing:
└── module9_combinedLR_stats_functions.R
    ├── calculate_summary_stats()
    ├── calculate_ratio_stats()
    ├── calculate_cutoffs()
    └── calculate_proportions_exceeding_cutoffs()
```

### Module Descriptions

#### Foundation
- **LR_kinship_utility_functions.R**: Central utilities, constants, and kinship LR calculations
  - Dependencies: None (loads data files directly)
  - Used by: All other modules

#### Simulation Modules
1. **module1_allele_simulator.R**: Simulates individual alleles from population frequencies
   - Dependencies: `LR_kinship_utility_functions.R`
   - Function: `simulate_allele(locus, population, allele_frequency_data)`

2. **module2_STR_profile_simulator.R**: Generates complete STR profiles (29 loci)
   - Dependencies: `module1_allele_simulator.R`, `LR_kinship_utility_functions.R`
   - Function: `simulate_str_profile(loci_list, population, allele_frequency_data)`

3. **module3_related_individual_simulator.R**: Creates genetically related pairs
   - Dependencies: `module1`, `module2`, `LR_kinship_utility_functions.R`
   - Functions: 
     - `simulate_related_individual()` - single related individual
     - `simulate_individual_pair()` - focal + related pair
     - `simulate_multiple_pairs()` - batch generation

#### Analysis Modules
4. **module4_single_locus_LR.R**: Calculates single-locus likelihood ratios
   - Dependencies: `LR_kinship_utility_functions.R`
   - Function: `calculate_single_locus_lr()`
   - Uses kinship coefficients and allele frequencies to compute LRs

5. **module5_combined_LR.R**: Combines LRs across loci sets
   - Dependencies: `LR_kinship_utility_functions.R`
   - Function: `calculate_combined_lr()`
   - Multiplies single-locus LRs for different marker panels

#### Extended Modules
6. **module6_single_combo_pair_generator.R**: Batch pair generation wrapper
7. **module7_single_pop_focal_generator.R**: Family structure simulator
8. **module8_unrelated_pool_generator.R**: Unrelated individual pool generator
9. **module9_combinedLR_stats_functions.R**: Statistical analysis functions

### Wrapper Scripts & Job Management

#### Simulation Pipeline
1. **sim_pairs.R**: R script for pair simulation
   - Called by: `sim_pairs.sh`
   - Dependencies: modules 1-3

2. **sim_pairs.sh**: SLURM array job for parallel simulation
   - Array: 1-300 (5 populations × 6 relationships × 10 chunks)
   - Output: `pairs_*.csv` files

#### LR Calculation Pipeline
3. **lr_wrapper.R**: R script for single-locus LR calculation
   - Called by: `lr_wrapper.sh`
   - Dependencies: module 4

4. **lr_wrapper.sh**: SLURM array job for LR calculation
   - Configured by: `lr_submission.sh`
   - Output: `LR/LR_*.csv` files

5. **lr_submission.sh**: Generates file list and submits LR jobs

#### Combined LR Pipeline
6. **combined_lr_wrapper.R**: R script for multi-locus LR
   - Called by: `combined_lr.sh`
   - Dependencies: module 5

7. **combined_lr.sh**: SLURM array job for combined LR
   - Configured by: `combined_lr_submission.sh`
   - Output: `combined_LR/combined_LR_*.csv` files

#### Analysis & Visualization
8. **analyze_lr_outputs.R**: Aggregates and analyzes all LR results
   - Called by: `analyze_lr_outputs.sh`
   - Dependencies: module 9
   - Creates summary statistics and prepares data for plotting

9. **plots_matched.R**: Visualizes results with correct assumptions
10. **plots_mismatched.R**: Visualizes population/relationship mismatch effects
11. **plots_proportion_exceeding_cutoffs.R**: Threshold analysis plots

## Workflow

### Step 1: Simulate Genotype Pairs
```bash
sbatch code/sim_pairs.sh
```
- Generates 240,000 pairs (10,000 per population-relationship combination)
- 5 populations: AfAm, Cauc, Hispanic, Asian, all
- 6 relationships: parent_child, full_siblings, half_siblings, cousins, second_cousins, unrelated
- Output: 240 files in `output/pairs_*.csv`

### Step 2: Calculate Single-Locus LRs
```bash
bash code/lr_submission.sh
```
- Tests each pair against multiple hypotheses:
  - 6 tested relationships
  - 5 tested populations
- Calculates LR for each of 29 STR loci
- Output: 240 files in `output/LR/`

### Step 3: Calculate Combined LRs (Optional)
```bash
bash code/combined_lr_submission.sh
```
- Combines single-locus LRs across loci sets:
  - core_13: Original CODIS markers
  - identifiler_15: Identifiler kit
  - expanded_20: Expanded CODIS
  - supplementary: Additional markers
  - autosomal_29: All autosomal STRs
- Output: 240 files in `output/combined_LR/`

### Step 4: Analyze Results
```bash
sbatch code/analyze_lr_outputs.sh [output_dir]
```
- Aggregates all LR results
- Calculates summary statistics
- Identifies population match/mismatch effects
- Computes LR ratios and cutoffs
- Output: `output/lr_analysis_YYYYMMDD/`

### Step 5: Generate Plots
```bash
# For matched population/relationship
sbatch code/plots_matched.sh output/lr_analysis_YYYYMMDD

# For mismatched scenarios
sbatch code/plots_mismatched.sh output/lr_analysis_YYYYMMDD

# For threshold analysis
sbatch code/plots_proportion_exceeding_cutoffs.sh output/lr_analysis_YYYYMMDD
```

## Output Files

### Simulation Output
- `pairs_{population}_{relationship}_n{count}_chunk{num}_{date}.csv`
  - Columns: batch_id, pair_id, population, known_relationship, locus, focal_A1, focal_A2, ind2_A1, ind2_A2

### LR Output
- `LR_{population}_{relationship}_n{count}_chunk{num}_{date}.csv`
  - Columns: batch_id, pair_id, population, known_relationship, locus, tested_relationship, tested_population, LR

### Combined LR Output
- `combined_LR_{population}_{relationship}_n{count}_chunk{num}_{date}.csv`
  - Columns: batch_id, pair_id, population, known_relationship, loci_set, tested_relationship, tested_population, combined_LR

### Analysis Output
- `combined_LR_all.rds`: Complete dataset
- `combined_LR_summary_stats.csv`: Summary statistics
- `combined_LR_ratio_summary.csv`: Population mismatch effects
- `combined_LR_cutoffs.csv`: Decision thresholds
- `combined_LR_exceeding_cutoffs.csv`: Performance metrics
