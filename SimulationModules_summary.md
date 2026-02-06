# PODFRIDGE Simulation Pipeline

**Last Updated**: January 2026  
**Current Dataset**: 1,000,000 simulated pairs (500,000 related + 500,000 unrelated)

## Overview

The PODFRIDGE (POPulation Data For Relationship Inference using DNA Genotype Evidence) pipeline simulates STR (Short Tandem Repeat) profiles for genetically related individual pairs and calculates likelihood ratios (LRs) for kinship analysis. These simulations are essential for understanding how different population assumptions and relationship hypotheses affect forensic DNA kinship testing accuracy.

### Current Simulation Scale

- **Related Pairs**: 20,000 pairs per relationship type × 5 relationship types = 100,000 per population
- **Unrelated Pairs**: 100,000 pairs per population
- **Total**: 500,000 related + 500,000 unrelated = 1,000,000 pairs
- **Populations**: 5 (AfAm, Cauc, Hispanic, Asian, all)
- **Relationships**: 6 types (parent_child, full_siblings, half_siblings, cousins, second_cousins, unrelated)

## Directory Structure

```
podfridge_simulations/
├── code/                    # All R scripts and SLURM job scripts
│   ├── module*.R           # Core simulation and analysis modules
│   ├── sim_pairs.R         # Pair simulation R wrapper
│   ├── sim_pairs.sh        # SLURM job for related pairs
│   ├── sim_pairs_unrelated.sh  # SLURM job for unrelated pairs
│   ├── lr_wrapper.R        # Single-locus LR calculation wrapper
│   ├── lr_wrapper.sh       # SLURM job template for LR calc
│   ├── lr_submission.sh    # LR job submission script
│   ├── combined_lr_wrapper.R    # Combined LR calculation wrapper
│   ├── combined_lr.sh      # SLURM job template for combined LR
│   ├── combined_lr_submission.sh # Combined LR job submission
│   ├── analyze_lr_outputs.R     # Analysis aggregation script
│   ├── analyze_lr_outputs.sh    # SLURM wrapper for analysis
│   ├── plots_matched.R          # Matched scenario plots
│   ├── plots_mismatched.R       # Mismatched scenario plots
│   ├── plots_proportion_exceeding_cutoffs.R  # Threshold analysis
│   ├── plots_proportion_exceeding_cutoffs.sh # SLURM wrapper
│   ├── plots_mismatched.sh      # SLURM wrapper for mismatched plots
│   ├── simulation_analysis.sh   # SLURM wrapper for R Markdown report
│   └── run_pipeline_bare.sh     # Manual step-by-step reference (copy/paste each section)
│
├── analysis/                # R Markdown analysis scripts and reports
│   ├── simulation_analysis.Rmd       # Statistical analysis R Markdown
│   └── simulation_analysis_*.html    # Generated HTML reports (dated)
│
├── data/                    # Reference data files
│   ├── df_allelefreq_combined.csv    # Allele frequencies for 5 populations
│   ├── kinship_coefficients.csv      # IBD sharing probabilities (k0, k1, k2)
│   └── core_CODIS_loci.csv           # Loci set definitions
│
├── logs/                    # SLURM job logs and error files
│   ├── sim_pairs_*.out/err
│   ├── unrel_pairs_*.out/err
│   ├── lr_calc_*.out/err
│   ├── combined_lr_*.out/err
│   └── analyze_lr_*.out/err
│
└── output/                  # Generated data and results
    ├── pairs_*.csv          # Simulated genotype pairs (1000 files)
    ├── LR/                  # Single-locus likelihood ratios (1000 files)
    │   └── LR_*.csv
    ├── combined_LR/         # Multi-locus combined LRs (1000 files)
    │   └── combined_LR_*.csv
    └── lr_analysis_*/       # Analysis results and plots
        ├── combined_LR_all.rds              # Complete combined dataset
        ├── combined_LR_match.csv            # Strictly matched data
        ├── combined_LR_mismatch.csv         # Mismatched data
        ├── combined_LR_summary_stats.csv    # Summary statistics
        ├── combined_LR_ratio_summary.csv    # Ratio statistics
        ├── combined_LR_ratios_raw.csv       # Raw ratio data
        ├── plots_matched/                   # Correct assumption plots
        ├── plots_mismatched/                # Incorrect assumption plots
        ├── plots_exceeding_cutoffs/         # Threshold analysis plots
        └── analysis_results/                # R Markdown report outputs
            ├── simulation_analysis_*.html   # Statistical report
            ├── classification_summary_29loci.png
            ├── proportions_with_classification.csv
            ├── detailed_rates_29loci.csv
            ├── stat_test_loci_effect.csv
            ├── stat_test_population_effect.csv
            └── stat_test_linear_trend.csv
```

## Module Architecture & Dependencies

### Core Simulation Modules (R Scripts)

```
Foundation:
└── LR_kinship_utility_functions.R
    ├── Provides: FALLBACK_FREQ constant (5/(2*1036) ≈ 0.002414)
    ├── Loads: df_allelefreq, loci_lists, kinship_matrix
    ├── Helper functions: 
    │   ├── count_shared_alleles() - counts IBD alleles with multiplicity
    │   ├── label_and_genotype() - deterministic genotype labeling
    │   ├── fetch_freq() - retrieves allele frequencies with fallback
    │   └── calculate_likelihood_ratio() - core LR calculation
    └── Main function: kinship_calculation() - long-format LR computation

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
  - Dependencies: Loads data files directly (`data/*.csv`)
  - Used by: All other modules
  - Key Constants:
    - `FALLBACK_FREQ = 5 / (2 * 1036)` for zero-frequency alleles
    - `loci_lists`: Named list of 5 loci sets (core_13, identifiler_15, expanded_20, supplementary, autosomal_29)
  - Key Functions:
    - `count_shared_alleles()`: Counts shared alleles WITH multiplicity (handles homozygous genotypes correctly)
    - `label_and_genotype()`: Creates deterministic allele-to-letter mapping (shared alleles → A, B, etc.)
    - `calculate_likelihood_ratio()`: Implements IBD-based LR formulas for 0, 1, or 2 shared alleles
    - `kinship_calculation()`: Main driver function returning long-format LR results

#### Simulation Modules
1. **module1_allele_simulator.R**: Simulates individual alleles from population frequencies
   - Dependencies: `LR_kinship_utility_functions.R`
   - Function: `simulate_allele(locus, population, allele_frequency_data)`
   - Features:
     - Applies FALLBACK_FREQ to zero-frequency alleles
     - Normalizes frequencies to sum to 1.0
     - Returns single allele as character string

2. **module2_STR_profile_simulator.R**: Generates complete STR profiles (29 loci)
   - Dependencies: `module1_allele_simulator.R`, `LR_kinship_utility_functions.R`
   - Function: `simulate_str_profile(loci_list, population, allele_frequency_data)`
   - Features:
     - Calls allele simulator twice per locus (independent draws)
     - Returns data frame with columns: population, locus, A1, A2

3. **module3_related_individual_simulator.R**: Creates genetically related pairs
   - Dependencies: `module1`, `module2`, `LR_kinship_utility_functions.R`
   - Functions: 
     - `calculate_shared_allele_probability()`: Samples IBD sharing pattern (none/one/both) based on k0/k1/k2
     - `simulate_related_individual()`: Generates related individual based on focal's genotype and relationship
     - `simulate_individual_pair()`: Creates focal + related pair in wide format
     - `simulate_multiple_pairs()`: Batch generation with unique pair IDs
   - Features:
     - Implements Mendelian inheritance patterns via IBD sharing
     - For homozygous focal (14-14), sampling doesn't matter as both alleles identical
     - Pair IDs formatted with leading zeros (001, 002, etc.)

#### Analysis Modules
4. **module4_single_locus_LR.R**: Calculates single-locus likelihood ratios
   - Dependencies: `LR_kinship_utility_functions.R`
   - Function: `calculate_single_locus_lr(pair_data, tested_relationship, tested_populations, ...)`
   - Features:
     - Uses kinship coefficients (k0, k1, k2) for tested relationship
     - Tests multiple population hypotheses per pair
     - Handles NA alleles gracefully (returns LR = NA)
     - Returns long-format results (one row per locus-population combination)

5. **module5_combined_LR.R**: Combines LRs across loci sets
   - Dependencies: `LR_kinship_utility_functions.R`
   - Function: `calculate_combined_lr(single_locus_results, loci_sets)`
   - Features:
     - Multiplies single-locus LRs within each loci set
     - Uses `prod(LR, na.rm = TRUE)` - silently handles NAs
     - Returns one combined LR per pair-loci_set-hypothesis combination
     - Groups by: batch_id, pair_id, population, known_relationship, tested_relationship, tested_population

#### Extended Modules
6. **module6_single_combo_pair_generator.R**: Batch pair generation wrapper
   - Functions:
     - `generate_pair_batch()`: Generates pairs for single pop-relationship combo, saves to CSV
     - `generate_multiple_pair_batches()`: Generates multiple combinations at once
     - `create_pair_combinations()`: Helper to create standard combinations
   - Note: Currently not used in main pipeline (replaced by direct sim_pairs.R calls)

7. **module7_single_pop_focal_generator.R**: Family structure simulator
   - Functions:
     - `generate_single_pop_focal()`: Creates focal individuals with flexible family structures
     - `generate_families_by_relationships()`: Internal function for family generation
     - `generate_multiple_pop_focal()`: Multi-population focal family generation
     - `create_standard_family_structure()`: Helper for common family types
   - Note: For future focal pool analysis (not currently implemented)

8. **module8_unrelated_pool_generator.R**: Unrelated individual pool generator
   - Functions:
     - `generate_unrelated_pool()`: Creates pool of unrelated individuals for single population
     - `generate_unrelated_individuals()`: Internal function using module2
     - `generate_multiple_pop_unrelated()`: Multi-population unrelated pool generation
   - Note: For future unrelated pool analysis (not currently implemented)

9. **module9_combinedLR_stats_functions.R**: Statistical analysis functions
   - Functions:
     - `calculate_summary_stats()`: Groups by all experimental factors, computes mean/median/sd/quantiles
     - `calculate_ratio_stats()`: Calculates LR_wrong/LR_correct ratios for population mismatch analysis
     - `calculate_cutoffs()`: Determines LR thresholds for fixed FPRs (1%, 0.1%, 0.01%) using unrelated pairs
     - `calculate_proportions_exceeding_cutoffs()`: Computes true positive rates at each cutoff
   - Features:
     - Filters for correct population (is_correct_pop == TRUE) when calculating cutoffs

### Wrapper Scripts & Job Management

#### Simulation Pipeline
1. **sim_pairs.R**: R script for pair simulation
   - Called by: `sim_pairs.sh` and `sim_pairs_unrelated.sh`
   - Dependencies: modules 1-3
   - Usage: `Rscript code/sim_pairs.R <POP> <RELATIONSHIP> <N_PAIRS> [CHUNK_NUM]`
   - Features:
     - Accepts population, relationship, n_pairs, and optional chunk number
     - Uses `simulate_multiple_pairs()` from module3
     - Incorporates chunk number into pair_id for global uniqueness
     - Outputs to `output/pairs_{pop}_{rel}_n{n}_chunk{chunk}_{date}.csv`

2. **sim_pairs.sh**: SLURM array job for parallel simulation of related pairs
   - Array: 1-600 (5 populations × 6 relationships × 20 chunks)
   - Default: 1000 pairs per chunk = 20,000 pairs per pop-relationship combo
   - Features:
     - Configurable N_PAIRS parameter (default 1000): `sbatch code/sim_pairs.sh [N_PAIRS]`
     - CHUNKS_PER_COMBO = 20 (creates chunks 1-20)
     - Calculates pop/relationship from array task ID
     - Each chunk generates unique pair IDs (c01_001, c01_002, ..., c20_1000)
   - Resources: 768 MB RAM, 30 min per task
   - Output: 300 files for pairs in `output/pairs_*.csv`

3. **sim_pairs_unrelated.sh**: SLURM array job for unrelated pair simulation
   - Array: 1-400 (5 populations × 80 chunks)
   - Default: 1000 pairs per chunk = 80,000 additional unrelated pairs per population
   - Features:
     - Configurable N_PAIRS parameter: `sbatch code/sim_pairs_unrelated.sh [N_PAIRS]`
     - CHUNKS_PER_POP = 80
     - Starts at chunk 21 (after related pairs chunks 1-20)
     - Creates chunks 21-100 for each population
   - Resources: 768 MB RAM, 2 hrs per task
   - Output: 400 files in `output/pairs_*.csv`
   - Total unrelated pairs: 100,000 per population (20k from sim_pairs.sh + 80k from this script)

#### LR Calculation Pipeline
4. **lr_wrapper.R**: R script for single-locus LR calculation
   - Called by: `lr_wrapper.sh`
   - Dependencies: module 4
   - Usage: `Rscript code/lr_wrapper.R <PAIRS_CSV_FILE>`
   - Features:
     - Tests all 6 relationship hypotheses (parent_child, full_siblings, half_siblings, cousins, second_cousins, unrelated)
     - Tests all 5 population hypotheses (AfAm, Cauc, Hispanic, Asian, all)
     - Outputs to `output/LR/LR_{filename}.csv`

5. **lr_wrapper.sh**: SLURM array job for LR calculation
   - Configured by: `lr_submission.sh` (sets array size dynamically)
   - Resources: 2 GB RAM, 4 hrs per task
   - Reads specific file from `output/lr_file_list.txt` based on array task ID
   - Output: `output/LR/LR_*.csv` files

6. **lr_submission.sh**: LR job submission script
   - Usage: `bash code/lr_submission.sh [CHUNK_RANGE]`
   - Features:
     - Generates file list of pairs to process
     - Supports chunk range filtering: 
       - Default (*): processes all chunks
       - `1..20`: processes chunks 1-20
       - `21..100`: processes chunks 21-100
     - Creates `output/lr_file_list.txt`
     - Updates lr_wrapper.sh with correct array size
     - Submits SLURM job and reports job ID
   - Chunk Range Examples:
     - `bash code/lr_submission.sh` - processes all chunks
     - `bash code/lr_submission.sh 1..20` - first 20k pairs per combo
     - `bash code/lr_submission.sh 21..100` - unrelated 80k

#### Combined LR Pipeline
7. **combined_lr_wrapper.R**: R script for multi-locus LR
   - Called by: `combined_lr.sh`
   - Dependencies: module 5
   - Usage: `Rscript code/combined_lr_wrapper.R <LR_CSV_FILE>`
   - Features:
     - Calculates combined LR for all 5 loci sets
     - Adds matching flags (is_correct_rel, is_correct_pop)
     - Outputs to `output/combined_LR/combined_LR_{filename}.csv`

8. **combined_lr.sh**: SLURM array job for combined LR
   - Configured by: `combined_lr_submission.sh`
   - Resources: 1 GB RAM, 1 hr per task
   - Reads specific file from `output/combined_lr_file_list.txt`
   - Output: `output/combined_LR/combined_LR_*.csv` files

9. **combined_lr_submission.sh**: Combined LR job submission script
   - Usage: `bash code/combined_lr_submission.sh [CHUNK_RANGE] [CUSTOM_FILE_LIST]`
   - Features:
     - Auto-generates file list from LR/ directory if no custom list provided
     - Supports chunk range filtering (same as lr_submission.sh)
     - Creates `output/combined_lr_file_list.txt`
     - Updates combined_lr.sh with correct array size and file list
   - Advanced Usage:
     - `bash code/combined_lr_submission.sh - my_list.txt` - use custom file list
     - `bash code/combined_lr_submission.sh none my_list.txt` - same as above

#### Analysis & Visualization
10. **analyze_lr_outputs.R**: Aggregates and analyzes all LR results
    - Called by: `analyze_lr_outputs.sh`
    - Dependencies: module 9
    - Usage: `Rscript code/analyze_lr_outputs.R [OUTPUT_DIR]`
    - Features:
      - Reads all combined_LR files
      - Combines into single dataset
      - Creates three key datasets:
        - `all_combined`: Complete dataset (saved as .rds)
        - `all_combined_pop_match`: Correct population only (for summary stats)
        - `all_combined_strict_match`: Both pop AND relationship correct (for LR distributions)
      - Applies module 9 functions:
        - Summary statistics (mean, median, sd, quantiles)
        - Ratio statistics (wrong_LR / correct_LR)
    - Output Files:
      - `combined_LR_all.rds` (compressed)
      - `combined_LR_match.csv` (strict match only)
      - `combined_LR_mismatch.csv`
      - `combined_LR_summary_stats.csv`
      - `combined_LR_ratio_summary.csv`
      - `combined_LR_ratios_raw.csv`

11. **analyze_lr_outputs.sh**: SLURM wrapper for analysis
    - Usage: `sbatch code/analyze_lr_outputs.sh [OUTPUT_DIR]`
    - Resources: 96 GB RAM, 2 hrs
    - Note: High memory requirement due to combining 1000+ files

12. **plots_matched.R**: Visualizes results with correct assumptions
    - Usage: `Rscript code/plots_matched.R <input_dir> [output_dir]`
    - Dependencies: tidyverse, data.table, scales
    - Input: `combined_LR_match.csv` (strict match data)
    - Features:
      - Filters for is_correct_pop == TRUE AND known_relationship == tested_relationship
      - Creates combined PDF: `all_matched_plots.pdf`
      - Individual plots:
        - `lr_distributions_boxplot_matched.png`: Log-scale boxplots by relationship
        - `mean_combined_lr_matched.png`: Line plots showing mean LR trends across loci sets
        - `proportions_exceeding_cutoffs_matched.png`: Bar plots for threshold performance
        - `heatmap_proportions_fixed_cutoff_matched.png`: Heatmap visualization
    - Population colors: AfAm=red, Asian=blue, Cauc=green, Hispanic=purple, all=orange

13. **plots_mismatched.R**: Visualizes population/relationship mismatch effects
    - Usage: `Rscript code/plots_mismatched.R <input_dir> [output_dir]`
    - Called by: `plots_mismatched.sh`
    - Dependencies: tidyverse, data.table, scales
    - Input: `combined_LR_all.rds` (full dataset)
    - Features:
      - Creates relationship mismatch boxplot showing all tested vs. known combinations
      - Creates population mismatch boxplots for parent-child and full-sibling hypotheses
      - Generates multi-page PDF: `mismatched_pops_all_relationships_LRboxplots.pdf`
      - Individual PNG files for each tested relationship × tested population combination
    - Resources (via SLURM wrapper): 42 GB RAM, 1 hr

14. **plots_proportion_exceeding_cutoffs.R**: Threshold analysis plots
    - Usage: `Rscript code/plots_proportion_exceeding_cutoffs.R <input_dir> <output_dir>`
    - Called by: `plots_proportion_exceeding_cutoffs.sh`
    - Dependencies: tidyverse, data.table, scales, shades
    - Input: `combined_LR_all.rds`
    - Features:
      - Calculates cutoffs from unrelated pairs (is_correct_pop == TRUE)
      - Generates proportions exceeding analysis for parent-child and full-sibling hypotheses
      - Includes statistical testing:
        - Chi-square tests for loci effect (does number of loci matter?)
        - Chi-square tests for population effect (does ancestry matter?)
        - Linear trend analysis (how does performance scale with loci?)
      - Creates classification summary plot (true positives vs. false positives)
      - Outputs detailed rate tables and key findings summary
    - Output Files:
      - `population_mismatch_proportions_analysis.pdf`
      - `classification_summary_29loci.png`
      - `combined_LR_cutoffs.csv`
      - `combined_LR_exceeding_cutoffs.csv`
      - `proportions_with_classification.csv`
      - `stat_test_loci_effect.csv`
      - `stat_test_population_effect.csv`
      - `stat_test_linear_trend.csv`
      - `detailed_rates_29loci.csv`
      - `key_findings_summary.txt`
    - Resources (via SLURM wrapper): 35 GB RAM, 1 hr

15. **plots_mismatched.sh**: SLURM wrapper for mismatched plots
    - Usage: `sbatch code/plots_mismatched.sh <input_dir> [output_subdir]`
    - Resources: 42 GB RAM, 1 hr

16. **plots_proportion_exceeding_cutoffs.sh**: SLURM wrapper for threshold plots
    - Usage: `sbatch code/plots_proportion_exceeding_cutoffs.sh <input_dir> <output_dir>`
    - Resources: 35 GB RAM, 1 hr

#### Statistical Report Generation
17. **simulation_analysis.Rmd**: R Markdown report for comprehensive statistical analysis
    - Location: `analysis/simulation_analysis.Rmd`
    - Dependencies: tidyverse, data.table, scales, knitr, kableExtra
    - Input: `combined_LR_all.rds`
    - Features:
      - Filters to parent-child and full-sibling tested relationships only
      - Calculates classification performance (true positives, related FP, unrelated FP)
      - Statistical hypothesis testing:
        - Chi-square tests for loci effect by population and relationship
        - Chi-square tests for population effect at 29 loci
        - Linear trend analysis showing slope/R²/significance
    - Generates: HTML report with interactive tables and plots

18. **simulation_analysis.sh**: SLURM submission script for report generation
    - Usage: `sbatch code/simulation_analysis.sh <input_dir> [output_dir]`
    - Resources: 21 GB RAM, 30 min
    - Features:
      - Automatically adds date to output filename
      - Renders R Markdown to HTML
      - Saves to `output/<output_dir>/analysis_results/`
    - Output: `simulation_analysis_YYYYMMDD.html`

#### Pipeline Management Script
19. **run_pipeline_bare.sh**: Manual step-by-step command reference
    - Usage: **Do NOT run as a script** - copy/paste each section manually
    - Features:
      - Lists all pipeline commands in order
      - Includes verification commands for each step
      - Minimal comments - just the essential commands
      - Designed for interactive execution
    - Best for: Learning the pipeline, troubleshooting individual steps, manual control

## Workflow

### Full Pipeline (Start to Finish)

#### Step 1: Simulate Genotype Pairs
```bash
# Simulate related pairs (20k per pop-relationship combo)
sbatch code/sim_pairs.sh

# Simulate additional unrelated pairs (80k per population)
sbatch code/sim_pairs_unrelated.sh
```
**Results**:
- 300 related pair files (5 pops × 6 relationships × 10 chunks × 1000 pairs = 300,000 related pairs)
- 400 unrelated pair files (5 pops × 80 chunks × 1000 pairs = 400,000 additional unrelated pairs)
- Total: 500,000 related + 100,000 unrelated (from sim_pairs) + 400,000 unrelated = 1,000,000 pairs
- Output: 1000 files in `output/pairs_*.csv`

**File Naming**:
- `pairs_{population}_{relationship}_n1000_chunk{01-20}_{YYYYMMDD}.csv` (related)
- `pairs_{population}_unrelated_n1000_chunk{21-100}_{YYYYMMDD}.csv` (unrelated)

**Computational Cost**:
- Time: ~3 minutes per chunk
- Memory: ~170 MB per task
- Total: ~50 CPU-hours for all simulations

#### Step 2: Calculate Single-Locus LRs
```bash
# Process all chunks (default)
bash code/lr_submission.sh

# Or process specific chunk ranges:
# bash code/lr_submission.sh 1..20     # Related pairs only
# bash code/lr_submission.sh 21..100  # Unrelated pairs only
```
**Process**:
- Tests each pair against 6 relationship hypotheses × 5 population hypotheses = 30 LR calculations per locus
- Each pair has 29 loci
- Per file: 1,000 pairs × 29 loci × 30 tests = 870,000 LR values

**Results**:
- Output: 1000 files in `output/LR/LR_*.csv`
- Each file: ~870,000 rows

**Computational Cost**:
- Time: ~1 hour per file (longest: 73 minutes)
- Memory: ~560 MB per task
- Total: ~1000 CPU-hours for all LR calculations

#### Step 3: Calculate Combined LRs
```bash
# Process all chunks (default)
bash code/combined_lr_submission.sh

# Or process specific chunk ranges:
# bash code/combined_lr_submission.sh 1..20
```
**Process**:
- Multiplies single-locus LRs across 5 loci sets
- Per file: 1,000 pairs × 5 loci sets × 6 tested relationships × 5 tested populations = 150,000 combined LR values

**Results**:
- Output: 1000 files in `output/combined_LR/combined_LR_*.csv`
- Each file: ~150,000 rows

**Computational Cost**:
- Time: ~5 seconds per file
- Memory: ~530 MB per task
- Total: ~1.5 CPU-hours for all combined LR calculations

#### Step 4: Analyze Results
```bash
sbatch code/analyze_lr_outputs.sh output/lr_analysis_YYYYMMDD
```
**Process**:
- Reads all 1000 combined_LR files
- Aggregates into single dataset (~150 million rows total)
- Filters and processes:
  - Strict match (pop AND relationship correct): for LR distribution plots
  - Population match only: for summary statistics
  - All data: for ratio analysis
- Applies statistical functions from module 9

**Results**:
- Output directory: `output/lr_analysis_YYYYMMDD/`
- Files created:
  - `combined_LR_all.rds` (~GB compressed)
  - `combined_LR_match.csv` (strict matches only)
  - `combined_LR_mismatch.csv` (mismatches only)
  - `combined_LR_summary_stats.csv`
  - `combined_LR_ratio_summary.csv`
  - `combined_LR_ratios_raw.csv`

**Computational Cost**:
- Time: ~20 minutes
- Memory: 96 GB (due to large dataset)
- Requires high-memory node

#### Step 5: Generate Plots
```bash
# Matched scenarios (correct population AND relationship)
Rscript code/plots_matched.R output/lr_analysis_YYYYMMDD lr_analysis_YYYYMMDD/plots_matched

# Mismatched scenarios (incorrect population or relationship)
sbatch code/plots_mismatched.sh output/lr_analysis_YYYYMMDD lr_analysis_YYYYMMDD/plots_mismatched

# Threshold analysis (proportion exceeding cutoffs)
sbatch code/plots_proportion_exceeding_cutoffs.sh lr_analysis_YYYYMMDD lr_analysis_YYYYMMDD/plots_exceeding_cutoffs
```

**Matched Plots** (correct assumptions):
- `all_matched_plots.pdf` (combined)
- `lr_distributions_boxplot_matched.png`
- `mean_combined_lr_matched.png`
- `proportions_exceeding_cutoffs_matched.png`
- `heatmap_proportions_fixed_cutoff_matched.png`

**Mismatched Plots** (incorrect assumptions):
- `mismatched_pops_all_relationships_LRboxplots.pdf`
- `relationship_mismatch_LRboxplot.png`
- Individual PNGs for each tested relationship × population

**Threshold Analysis Plots**:
- `population_mismatch_proportions_analysis.pdf`
- `classification_summary_29loci.png`
- Statistical test CSV files
- Key findings summary

**Computational Cost**:
- plots_matched: ~5 minutes, ~10 GB RAM
- plots_mismatched: ~5 minutes, 42 GB RAM
- plots_proportion_exceeding_cutoffs: ~2 minutes, 35 GB RAM

#### Step 6: Generate Statistical Report
```bash
sbatch code/simulation_analysis.sh lr_analysis_YYYYMMDD
```
**Process**:
- Renders R Markdown report with comprehensive statistical analysis
- Performs chi-square tests and trend analysis
- Generates classification tables and visualizations

**Results**:
- `simulation_analysis_YYYYMMDD.html` (interactive report)
- CSV files with statistical test results
- Classification performance plots

**Computational Cost**:
- Time: ~7 minutes
- Memory: ~17 GB RAM

### Alternative Workflows

#### Quick Reference:
**For Manual Step-by-Step Execution:**
```bash
# DO NOT RUN AS: bash code/run_pipeline_bare.sh
# Instead, open the file and copy/paste each section one at a time:
cat code/run_pipeline_bare.sh

# Copy and paste each section, wait for completion before next section
# Example:
sbatch code/sim_pairs.sh
sbatch code/sim_pairs_unrelated.sh
# Wait: squeue -u $USER shows no jobs
# Verify: ls output/pairs_*.csv | wc -l
# Then proceed to next section...
```

#### Testing with Small Dataset
```bash
# Simulate small test dataset (100 pairs per chunk, 3 chunks)
sbatch --array=1-3 code/sim_pairs.sh 100

# Process just those chunks
bash code/lr_submission.sh 1..3
bash code/combined_lr_submission.sh 1..3

# Analyze
sbatch code/analyze_lr_outputs.sh output/lr_analysis_test
```

#### Incremental Updates
```bash
# Add more unrelated pairs (additional 500 per chunk)
sbatch code/sim_pairs_unrelated.sh 500

# Process only the new chunks
bash code/lr_submission.sh 101..150
bash code/combined_lr_submission.sh 101..150

# Re-run analysis with all data
sbatch code/analyze_lr_outputs.sh output/lr_analysis_updated
```

## Output Files

### Simulation Output
**Filename Pattern**: `pairs_{population}_{relationship}_n{count}_chunk{num}_{date}.csv`

**Example**: `pairs_AfAm_full_siblings_n1000_chunk05_20251219.csv`

**Columns**:
- `batch_id`: Timestamp when chunk was created (YYYYMMDD_HHMMSS)
- `pair_id`: Unique pair identifier incorporating chunk (e.g., "c05_042")
- `population`: True population (AfAm, Cauc, Hispanic, Asian, all)
- `known_relationship`: True relationship between individuals
- `locus`: STR marker name (29 loci total)
- `focal_A1`, `focal_A2`: Focal individual's alleles
- `ind2_A1`, `ind2_A2`: Related individual's alleles

**Size**: ~1.8 MB per 1000-pair file (29,001 rows including header)

### LR Output
**Filename Pattern**: `LR_{population}_{relationship}_n{count}_chunk{num}_{date}.csv`

**Example**: `LR_AfAm_full_siblings_n1000_chunk05_20251219.csv`

**Columns**:
- All columns from pairs file, plus:
- `shared_alleles`: Number of IBD alleles (0, 1, or 2)
- `genotype_match`: Genotype pattern (e.g., "AA-AB", "AB-AB")
- `tested_relationship`: Hypothesis being tested
- `tested_population`: Allele frequencies used
- `LR`: Single-locus likelihood ratio

**Size**: ~50 MB per file (870,000 rows)

### Combined LR Output
**Filename Pattern**: `combined_LR_{population}_{relationship}_n{count}_chunk{num}_{date}.csv`

**Example**: `combined_LR_AfAm_full_siblings_n1000_chunk05_20251219.csv`

**Columns**:
- `batch_id`, `pair_id`: Pair identifiers
- `population`: True population
- `known_relationship`: True relationship
- `loci_set`: Marker panel used (core_13, identifiler_15, expanded_20, supplementary, autosomal_29)
- `tested_relationship`: Hypothesis being tested
- `tested_population`: Allele frequencies used
- `combined_LR`: Product of single-locus LRs
- `is_correct_rel`: Boolean, known == tested relationship
- `is_correct_pop`: Boolean, population == tested_population

**Size**: ~8 MB per file (150,000 rows)

### Analysis Output

**Primary Dataset**:
- `combined_LR_all.rds`: Complete aggregated dataset (compressed, ~GB size)
  - All 1000 files combined
  - ~150 million rows total

**Filtered Datasets**:
- `combined_LR_match.csv`: Strict match (pop AND relationship correct)
  - Used for: LR distribution plots showing expected performance
  - Filter: `is_correct_pop == TRUE & known_relationship == tested_relationship`
  
- `combined_LR_mismatch.csv`: Any mismatch
  - Used for: Ratio analysis and mismatch effect plots
  - Filter: `is_correct_pop == FALSE`

**Summary Statistics**:
- `combined_LR_summary_stats.csv`
  - Groups by: known_relationship, population, loci_set, tested_population, tested_relationship, is_correct_pop
  - Metrics: n, mean_LR, median_LR, sd_LR, min_LR, max_LR, lower_95, upper_95

**Ratio Analysis**:
- `combined_LR_ratio_summary.csv`
  - Ratio = LR_wrong_pop / LR_correct_pop
  - Groups by: population, known_relationship, tested_relationship, loci_set, tested_population
  - Metrics: n, mean_ratio, median_ratio, sd_ratio, min_ratio, max_ratio, lower_95, upper_95

- `combined_LR_ratios_raw.csv`
  - Individual pair-level ratios (not summarized)
  - Used for detailed analysis and visualization

### Statistical Report Output
Generated in `output/lr_analysis_YYYYMMDD/analysis_results/`:

**Main Report**:
- `simulation_analysis_YYYYMMDD.html`: Interactive HTML report with:
  - Data summaries and distribution plots
  - Classification performance analysis
  - Statistical test results with interpretation
  - Key findings and conclusions

**Supporting Files**:
- `classification_summary_29loci.png`: Visual summary of classification performance
- `proportions_with_classification.csv`: Full proportions data with TP/FP labels
- `detailed_rates_29loci.csv`: Detailed rates table at maximum loci
- `stat_test_loci_effect.csv`: Chi-square test results for loci effect
- `stat_test_population_effect.csv`: Chi-square test results for population effect  
- `stat_test_linear_trend.csv`: Linear regression results showing trends with loci


## Key Technical Details

### Unique Pair Identification
- **Global Uniqueness**: Achieved via `(batch_id, pair_id)` combination
  - `batch_id`: Timestamp when chunk was created (unique per simulation run)
  - `pair_id`: Incorporates chunk number (e.g., "c05_042" = chunk 5, pair 42)

### Memory Management
- **High-memory steps**:
  - analyze_lr_outputs.R: 96 GB (combines 150M rows)
  - plots_mismatched: 42 GB
  - plots_proportion_exceeding_cutoffs: 35 GB

### Chunk Range Syntax
- **Correct**: `1..20` (double dots, compatible with bash brace expansion)
- **Incorrect**: `1-20` (single dash causes error)
- **Validation**: Scripts check for incorrect single dash and provide helpful error message

## Computational Resources Summary

**Total for 1M Pair Pipeline**:
- **Simulation**: ~3 minutes, 155 MB RAM per task
- **LR Calculation**: ~1.5 hours, ~400 MB RAM per task
- **Combined LR**: ~1 minute, ~95 MB RAM per task
- **Analysis**: ~20 minutes, 96 GB RAM (single task)
- **Plotting**: ~20 minutes total, 25-36 GB RAM
- **Statistical Report**: ~5 minutes, 10 GB RAM

**Storage Requirements**:
- Simulation files: ~1 GB (1000 files × 1 MB)
- LR files: ~16 GB (1000 files × 16 MB)
- Combined LR files: ~3.3 GB (1000 files × 3.3 MB)
- Analysis outputs: ~7.6 GB
- Total: ~27.3 GB for complete pipeline

**Timeline** (with parallelization):
- Step 1 (Simulation): ~3 hour wall time (array parallelization)
- Step 2 (LR Calc): ~1.5 hours wall time (array parallelization)
- Step 3 (Combined LR): ~1 minutes wall time (array parallelization)
- Step 4 (Analysis): ~20 minutes
- Step 5 (Plotting): ~15-20 minutes
- Step 6 (Report): ~10 minutes
- **Total**: ~5 hours from start to final report

---

**Last Updated**: January 29, 2026
