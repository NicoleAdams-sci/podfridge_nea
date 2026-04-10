# PODFRIDGE Simulation Pipeline

**Last Updated**: April 2026  
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
│   ├── prepare_combined_lr_intermediates.R  # Intermediate CSV prep (NEW)
│   ├── prepare_combined_lr_intermediates.sh # SLURM wrapper for above (NEW)
│   ├── plots_matched_publication.R          # Publication matched plots (NEW)
│   ├── plots_mismatched_population.R        # Population mismatch plots (NEW)
│   ├── plots_mismatched_relationship.R      # Relationship discrimination plots (NEW)
│   ├── plots_cutoffs_publication.R          # Classification/FPR threshold plots (NEW)
│   ├── run_statistical_tests.R              # All inferential tests (NEW)
│   ├── run_statistical_tests.sh             # SLURM wrapper for stats (NEW)
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
│   ├── analyze_lr_*.out/err
│   ├── prep_intermediates_*.out/err
│   └── stats_tests_*.out/err
│
└── output/                  # Generated data and results
    ├── pairs_*.csv          # Simulated genotype pairs (1000 files)
    ├── LR/                  # Single-locus likelihood ratios (1000 files)
    │   └── LR_*.csv
    ├── combined_LR/         # Multi-locus combined LRs (1000 files)
    │   └── combined_LR_*.csv
    └── lr_analysis_*/       # Analysis results and plots
        ├── combined_LR_all.rds              # Complete combined dataset
        ├── combined_LR_match.csv.gz         # Strictly matched data (gzip compressed)
        ├── combined_LR_mismatch.csv         # Mismatched data
        ├── combined_LR_summary_stats.csv    # Summary statistics
        ├── combined_LR_ratio_summary.csv    # Ratio statistics
        ├── combined_LR_ratios_raw.csv       # Raw ratio data
        ├── proportions_with_classification.csv  # Intermediate: cutoff/FPR data
        ├── mismatched_pop_robustness.csv        # Intermediate: pop mismatch lines
        ├── mismatched_pop_heatmap.csv           # Intermediate: pop mismatch heatmap
        ├── plots_matched/                   # Publication matched plots
        ├── plots_mismatched_population/     # Population mismatch plots
        ├── plots_mismatched_relationship/   # Relationship discrimination plots
        ├── plots_cutoffs/                   # Classification/FPR threshold plots
        ├── stats/                           # Statistical test results (NEW)
        │   ├── stats_matched.csv
        │   ├── stats_mismatched_population.csv
        │   └── stats_fpr_cutoffs.csv
        └── analysis_results/                # R Markdown report outputs
            ├── simulation_analysis_*.html
            └── ...
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

#### Analysis & Intermediate Preparation

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
      - `combined_LR_match.csv.gz` (gzip compressed strict match data)
      - `combined_LR_mismatch.csv`
      - `combined_LR_summary_stats.csv`
      - `combined_LR_ratio_summary.csv`
      - `combined_LR_ratios_raw.csv`

11. **analyze_lr_outputs.sh**: SLURM wrapper for analysis
    - Usage: `sbatch code/analyze_lr_outputs.sh [OUTPUT_DIR]`
    - Resources: 96 GB RAM, 2 hrs
    - Note: High memory requirement due to combining 1000+ files

12. **prepare_combined_lr_intermediates.R** *(NEW)*: Generates aggregated intermediate CSVs for plotting scripts
    - Called by: `prepare_combined_lr_intermediates.sh`
    - Usage: `Rscript code/prepare_combined_lr_intermediates.R <input_dir> [output_dir]`
    - Purpose: Loads `combined_LR_all.rds` once and writes three small aggregated CSVs, avoiding
      the need for each plotting script to reload the full dataset independently
    - Features:
      - Section 1 — Pre-filter diagnostics: reports zero LR counts, NA rates, exclusion rates by
        population match status; output goes to SLURM log for post-hoc inspection
      - Section 2 — Preprocessing: adds `pop_match`, `pop_match_label`, `rel_match`, `log10_LR` columns
      - Section 3 — `proportions_with_classification.csv`: empirical FPR cutoffs from unrelated pairs
        (1%, 0.1%, 0.01%), proportions exceeding each cutoff, and classification labels
        (True Positive / Related FP / Unrelated FP); convenience columns `prop_LR_gt_1/10/100/1000`
      - Section 4 — `mismatched_pop_robustness.csv`: true positive pairs (rel_match == TRUE, parent_child
        and full_siblings only); median and mean log10(LR) per population × tested_population ×
        relationship × loci_set × pop_match_label cell
      - Section 5 — `mismatched_pop_heatmap.csv`: pair-level log10(LR_wrong/LR_correct) joined by
        pair ID; diagonal rows added with log10_ratio = 0 by construction; `is_diagonal` flag retained
    - Output Files (written to `<output_dir>/`):
      - `proportions_with_classification.csv`
      - `mismatched_pop_robustness.csv`
      - `mismatched_pop_heatmap.csv`

13. **prepare_combined_lr_intermediates.sh** *(NEW)*: SLURM wrapper for intermediate preparation
    - Usage: `sbatch code/prepare_combined_lr_intermediates.sh <input_dir>`
    - Resources: 48 GB RAM, 1 hr
    - Validates `combined_LR_all.rds` exists before running
    - Prints next-step commands on success

#### Publication Plotting Scripts *(ALL NEW / REDESIGNED)*

> **Note**: The previous plotting scripts (`plots_matched.R`, `plots_mismatched.R`,
> `plots_proportion_exceeding_cutoffs.R`) have been replaced with publication-quality
> scripts designed for manuscript figures. The new scripts use the
> **Okabe-Ito colorblind-safe palette** throughout and are typically run interactively
> (`Rscript`) rather than via SLURM.

**Shared color palettes (Okabe-Ito, consistent across all publication scripts)**:
- Relationship colors: Parent-Child=#D55E00, Full Siblings=#E69F00, Half Siblings=#56B4E9,
  Cousins=#009E73, Second Cousins=#CC79A7, Unrelated=#999999
- Population colors: AfAm=#0072B2, Asian=#009E73, Cauc=#56B4E9, Hispanic=#CC79A7, all=#999999

14. **plots_matched_publication.R** *(replaces plots_matched.R)*: Matched scenario publication figures
    - Usage: `Rscript code/plots_matched_publication.R <input_dir> [output_dir]`
    - Input: `combined_LR_match.csv.gz` (reads directly; does not require intermediates step)
    - Dependencies: tidyverse, data.table, scales, patchwork
    - Figures:
      - **Main text figure** (`matched_main_lr_distributions.pdf/.png`): Violin plots of log10(LR)
        by relationship type across all 5 loci panels; `population = "all"` only; quantile lines +
        mean diamonds; faceted by loci set (5 columns); 12 × 5 inches
      - **Supplement figure** (`matched_supp_lr_by_population.pdf/.png`): Grouped boxplots showing
        all 5 populations per relationship type; faceted by loci set (3 columns); 12 × 10 inches
    - Descriptive Output: `matched_summary_statistics_for_ms.csv` (mean, median, SD, IQR per
      relationship × loci set for the "all" population)
    - Prints console summaries for parent-child and unrelated LRs as a quick sanity check

15. **plots_mismatched_population.R** *(new — replaces part of plots_mismatched.R)*: Population mismatch figures
    - Usage: `Rscript code/plots_mismatched_population.R <input_dir> [output_dir]`
    - Input: `mismatched_pop_robustness.csv` + `mismatched_pop_heatmap.csv` (from intermediates step)
    - Dependencies: tidyverse, data.table, scales, patchwork
    - Figures:
      - **Main figure** (`mismatched_population_robustness.png`): Line plots of median log10(LR) vs.
        number of loci for matched vs. mismatched population frequencies (solid/dashed lines);
        faceted by true population (rows) × relationship type (columns); 10 × 8 inches
      - **Supplement figure** (`mismatched_population_supp_heatmap.png`): Tile heatmap of
        median log10(LR_wrong / LR_correct) for each true × tested population combination;
        gold border highlights diagonal (no mismatch); color scale capped at 3 (10³× inflation);
        faceted by loci set (rows) × relationship (columns); 12 × 20 inches
    - Summary output: `mismatched_population_comparison_detailed.csv`

16. **plots_mismatched_relationship.R** *(new — replaces part of plots_mismatched.R)*: Relationship discrimination figure
    - Usage: `Rscript code/plots_mismatched_relationship.R <input_dir> [output_dir]`
    - Input: `combined_LR_all.rds` (loads the full dataset directly)
    - Dependencies: tidyverse, data.table, scales, patchwork
    - Figures:
      - **Main figure** (`mismatched_relationship_discrimination.png`): Boxplots of combined LR by
        tested relationship hypothesis for true pairs; parent-child and full-sibling hypotheses only;
        Core 13, Expanded 20, and Autosomal 29 loci panels only; gold shading on panels where
        tested hypothesis matches true relationship; colored by true population; 12 × 10 inches
    - Note: This script uses the legacy population color scheme (red/blue/green/purple/orange)
      rather than the Okabe-Ito scheme; retained because it visualises continuous LR
      distributions rather than threshold-based summaries

17. **plots_cutoffs_publication.R** *(replaces plots_proportion_exceeding_cutoffs.R)*: Classification performance figures
    - Usage: `Rscript code/plots_cutoffs_publication.R <input_dir> [output_dir]`
    - Input: `proportions_with_classification.csv` (from intermediates step)
    - Dependencies: tidyverse, data.table, scales, patchwork
    - Figures:
      - **Figure 1 — Main text** (`cutoff_classification_0.01fpr_29loci.png`): Grouped bar chart of
        proportion of pairs exceeding the 0.01% FPR threshold at 29 autosomal loci only;
        faceted by tested hypothesis (rows) × true relationship (columns); colored by population;
        18 × 10 inches
      - **Figure 2 — Supplement** (`cutoff_classification_fpr_29loci.png`): Same layout extended
        across all three FPR thresholds (1%, 0.1%, 0.01%); 6 row facets (2 hypotheses × 3 thresholds);
        classification (TP / Related FP / Unrelated FP) on x-axis; 18 × 16 inches
      - **Figure 3 — Supplement** (`cutoff_supp_heatmap_fp_rates.png`): Heatmap of false positive
        rates across all loci panels and all 4 thresholds (including Fixed LR > 1); stratified by
        population and tested hypothesis; two-panel layout (A: parent-child, B: full siblings);
        YlOrRd fill scale; 22 × 20 inches

#### Statistical Tests *(NEW — separated from plotting scripts)*

18. **run_statistical_tests.R** *(NEW)*: All inferential tests for manuscript results section
    - Usage: `Rscript code/run_statistical_tests.R <input_dir> [output_dir]`
    - Default output dir: `<input_dir>/stats/`
    - Dependencies: tidyverse, data.table, rstatix (kruskal_effsize)
    - Note: With large simulation N, p-values will be near-zero throughout; effect sizes are
      the informative statistic in all tests
    - **Section 1 — Matched population** (input: `combined_LR_match.csv.gz`):
      - Test 1: Kruskal-Wallis + epsilon-squared — does log10(LR) differ by relationship type,
        per loci set? (`population = "all"` only)
      - Test 2: Spearman correlation — does log10(LR) increase with loci count,
        per relationship type? (computed on per-group medians)
      - Test 3: Kruskal-Wallis + epsilon-squared — does population affect log10(LR) at
        Autosomal 29? (named populations only, per relationship)
      - Output: `stats_matched.csv`
    - **Section 2 — Mismatched population** (input: `combined_LR_all.rds`; filtered to true positives):
      - Test 4: Paired Wilcoxon signed-rank — does using wrong population frequencies inflate
        log10(LR)? Per true_population × tested_population × relationship × loci_set cell;
        effect size = median log10(LR_wrong / LR_correct); +1 = 10× inflation
      - Test 5: Linear model (ANOVA) — which factors drive LR inflation? Fitted on per-cell
        medians to avoid pseudoreplication; reports eta-squared per term
        (true_population, tested_population, loci_set) per relationship
      - Output: `stats_mismatched_population.csv`
    - **Section 3 — FPR cutoff analysis** (input: `proportions_with_classification.csv`):
      - Test 6: Cochran-Armitage trend test — does FPR change monotonically with loci count?
        Per tested_relationship × population × classification × threshold combination
      - Test 7: Logistic regression with drop1 LRT — do relationship type, population, and
        loci set affect FPR? Per tested_relationship × threshold; effect size = deviance explained
        (term deviance / null deviance)
      - Output: `stats_fpr_cutoffs.csv`
    - All three output CSVs use a unified long format with columns:
      section, test, question, grouping, statistic, df, p_value, effect_size, effect_size_type,
      magnitude, n, note

19. **run_statistical_tests.sh** *(NEW)*: SLURM wrapper for statistical tests
    - Usage: `sbatch code/run_statistical_tests.sh <input_dir>`
    - Resources: 64 GB RAM, 2 hrs
    - Validates all three required input files before running:
      `combined_LR_all.rds`, `combined_LR_match.csv.gz`, `proportions_with_classification.csv`

#### Statistical Report Generation
20. **simulation_analysis.Rmd**: R Markdown report for descriptive statistical analysis
    - Location: `analysis/simulation_analysis.Rmd`
    - Dependencies: tidyverse, data.table, scales, knitr, kableExtra
    - Input: `combined_LR_all.rds`
    - Note: Inferential tests have been moved to `run_statistical_tests.R`; this script
      now focuses on descriptive summaries, classification tables, and visualizations
    - Generates: HTML report with interactive tables and plots

21. **simulation_analysis.sh**: SLURM submission script for report generation
    - Usage: `sbatch code/simulation_analysis.sh <input_dir> [output_dir]`
    - Resources: 21 GB RAM, 30 min
    - Output: `simulation_analysis_YYYYMMDD.html`

#### Pipeline Management Script
22. **run_pipeline_bare.sh**: Manual step-by-step command reference
    - Usage: **Do NOT run as a script** - copy/paste each section manually
    - Includes verification commands for each step
    - Designed for interactive execution

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

**Computational Cost**:
- Time: ~3 minutes per chunk
- Memory: ~170 MB per task

#### Step 2: Calculate Single-Locus LRs
```bash
bash code/lr_submission.sh
```
**Process**: Tests each pair against 6 relationship hypotheses × 5 population hypotheses = 30 LR calculations per locus; 29 loci per pair = 870,000 LR values per file.

**Computational Cost**:
- Time: ~1 hour per file
- Memory: ~560 MB per task

#### Step 3: Calculate Combined LRs
```bash
bash code/combined_lr_submission.sh
```
**Process**: Multiplies single-locus LRs across 5 loci sets; 150,000 combined LR values per file.

**Computational Cost**:
- Time: ~5 seconds per file
- Memory: ~530 MB per task

#### Step 4: Analyze Results
```bash
sbatch code/analyze_lr_outputs.sh output/lr_analysis_YYYYMMDD
```
**Process**: Reads all 1000 combined_LR files; aggregates into single dataset (~150 million rows); 
filters to strict match and mismatch subsets; applies module 9 statistical functions.

**Output Files**:
- `combined_LR_all.rds` (~GB compressed)
- `combined_LR_match.csv.gz` (gzip compressed strict matches)
- `combined_LR_mismatch.csv`
- `combined_LR_summary_stats.csv`
- `combined_LR_ratio_summary.csv`
- `combined_LR_ratios_raw.csv`

**Computational Cost**: ~20 minutes, 96 GB RAM (high-memory node required)

#### Step 4.5: Prepare Intermediate Files *(NEW STEP)*
```bash
sbatch code/prepare_combined_lr_intermediates.sh output/lr_analysis_YYYYMMDD
```
**Purpose**: Loads `combined_LR_all.rds` once and pre-aggregates three small CSVs consumed by the
downstream publication plotting scripts. Running this step avoids each plotting script having to
independently reload the full dataset.

**Output Files** (written alongside other analysis outputs):
- `proportions_with_classification.csv` → consumed by `plots_cutoffs_publication.R`
- `mismatched_pop_robustness.csv` → consumed by `plots_mismatched_population.R`
- `mismatched_pop_heatmap.csv` → consumed by `plots_mismatched_population.R`

**Computational Cost**: ~10 minutes, 48 GB RAM

#### Step 5: Generate Publication Plots
```bash
# Matched scenarios — run interactively
Rscript code/plots_matched_publication.R output/lr_analysis_YYYYMMDD \
    output/lr_analysis_YYYYMMDD/plots_matched

# Population mismatch — run interactively
Rscript code/plots_mismatched_population.R output/lr_analysis_YYYYMMDD \
    output/lr_analysis_YYYYMMDD/plots_mismatched_population

# Relationship discrimination — run interactively (loads full .rds)
Rscript code/plots_mismatched_relationship.R output/lr_analysis_YYYYMMDD \
    output/lr_analysis_YYYYMMDD/plots_mismatched_relationship

# FPR/classification threshold plots — run interactively
Rscript code/plots_cutoffs_publication.R output/lr_analysis_YYYYMMDD \
    output/lr_analysis_YYYYMMDD/plots_cutoffs
```

**Matched Plots** (`plots_matched/`):
- `matched_main_lr_distributions.pdf/.png` (violin, all populations combined)
- `matched_supp_lr_by_population.pdf/.png` (boxplot, per population)
- `matched_summary_statistics_for_ms.csv`

**Population Mismatch Plots** (`plots_mismatched_population/`):
- `mismatched_population_robustness.png` (line plots, matched vs. mismatched freq)
- `mismatched_population_supp_heatmap.png` (heatmap of LR inflation)
- `mismatched_population_comparison_detailed.csv`

**Relationship Discrimination Plots** (`plots_mismatched_relationship/`):
- `mismatched_relationship_discrimination.png` (boxplots across hypothesis tests)

**Classification/FPR Plots** (`plots_cutoffs/`):
- `cutoff_classification_0.01fpr_29loci.png` (main text, single threshold)
- `cutoff_classification_fpr_29loci.png` (supplement, all thresholds)
- `cutoff_supp_heatmap_fp_rates.png` (supplement, FPR heatmap)

#### Step 6: Run Statistical Tests *(NEW STEP)*
```bash
sbatch code/run_statistical_tests.sh output/lr_analysis_YYYYMMDD
```
**Process**: Runs 7 inferential tests across 3 sections; writes unified long-format result tables.
Effect sizes are the primary reported statistic given the large simulation N.

**Output** (`output/lr_analysis_YYYYMMDD/stats/`):
- `stats_matched.csv` (Tests 1–3: KW, Spearman)
- `stats_mismatched_population.csv` (Tests 4–5: Wilcoxon, linear model)
- `stats_fpr_cutoffs.csv` (Tests 6–7: Cochran-Armitage, logistic regression)

**Computational Cost**: ~30–60 minutes, 64 GB RAM (Section 2 loads full .rds)

#### Step 7: Generate Statistical Report (Optional)
```bash
sbatch code/simulation_analysis.sh lr_analysis_YYYYMMDD
```
**Output**: `simulation_analysis_YYYYMMDD.html` — interactive descriptive report

**Computational Cost**: ~7 minutes, ~17 GB RAM

### Alternative Workflows

#### Quick Reference:
**For Manual Step-by-Step Execution:**
```bash
# DO NOT RUN AS: bash code/run_pipeline_bare.sh
# Instead, open the file and copy/paste each section one at a time:
cat code/run_pipeline_bare.sh
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

# Prepare intermediates
sbatch code/prepare_combined_lr_intermediates.sh output/lr_analysis_test

# Plots and stats
Rscript code/plots_matched_publication.R output/lr_analysis_test
sbatch code/run_statistical_tests.sh output/lr_analysis_test
```

## Output Files

### Simulation Output
**Filename Pattern**: `pairs_{population}_{relationship}_n{count}_chunk{num}_{date}.csv`

**Columns**:
- `batch_id`: Timestamp when chunk was created (YYYYMMDD_HHMMSS)
- `pair_id`: Unique pair identifier incorporating chunk (e.g., "c05_042")
- `population`: True population (AfAm, Cauc, Hispanic, Asian, all)
- `known_relationship`: True relationship between individuals
- `locus`: STR marker name (29 loci total)
- `focal_A1`, `focal_A2`: Focal individual's alleles
- `ind2_A1`, `ind2_A2`: Related individual's alleles

**Size**: ~1.8 MB per 1000-pair file

### Combined LR Output
**Columns**:
- `batch_id`, `pair_id`: Pair identifiers
- `population`: True population
- `known_relationship`: True relationship
- `loci_set`: Marker panel used
- `tested_relationship`: Hypothesis being tested
- `tested_population`: Allele frequencies used
- `combined_LR`: Product of single-locus LRs
- `is_correct_rel`: Boolean, known == tested relationship
- `is_correct_pop`: Boolean, population == tested_population

### Analysis Output

**Primary Dataset**:
- `combined_LR_all.rds`: Complete aggregated dataset (~150 million rows, compressed)

**Filtered Datasets**:
- `combined_LR_match.csv.gz`: Strict match (pop AND relationship correct) — **now gzip compressed**
- `combined_LR_mismatch.csv`: Any population mismatch

**Intermediate Aggregated Files** *(new)*:
- `proportions_with_classification.csv`: Proportions of pairs exceeding each LR threshold, with
  classification labels (True Positive / Related FP / Unrelated FP) and `n_loci` column
- `mismatched_pop_robustness.csv`: Median/mean log10(LR) per population × relationship × loci cell
  for matched and mismatched frequency usage
- `mismatched_pop_heatmap.csv`: Median log10(LR_wrong/LR_correct) per true × tested population cell;
  includes `is_diagonal` flag

**Statistical Test Results** *(new, in `stats/` subdirectory)*:
- `stats_matched.csv`: KW + Spearman tests (Tests 1–3)
- `stats_mismatched_population.csv`: Wilcoxon + linear model tests (Tests 4–5)
- `stats_fpr_cutoffs.csv`: Cochran-Armitage + logistic regression tests (Tests 6–7)
- Unified long format: section, test, question, grouping, statistic, df, p_value, effect_size,
  effect_size_type, magnitude, n, note

## Key Technical Details

### Unique Pair Identification
- **Global Uniqueness**: Achieved via `(batch_id, pair_id)` combination
  - `batch_id`: Timestamp when chunk was created (unique per simulation run)
  - `pair_id`: Incorporates chunk number (e.g., "c05_042" = chunk 5, pair 42)
- Required for correct pair-level joining in population mismatch analyses (Steps 4.5, 6)

### Memory Management
- **High-memory steps**:
  - `analyze_lr_outputs.R`: 96 GB (combines 150M rows)
  - `prepare_combined_lr_intermediates.R`: 48 GB
  - `run_statistical_tests.R`: 64 GB (loads full .rds for Section 2)
  - `plots_mismatched_relationship.R`: ~42 GB (loads full .rds)

### Chunk Range Syntax
- **Correct**: `1..20` (double dots, compatible with bash brace expansion)
- **Incorrect**: `1-20` (single dash causes error)
- Scripts validate for incorrect single dash and provide a helpful error message

## Computational Resources Summary

**Total for 1M Pair Pipeline**:
- **Simulation**: ~3 minutes, 155 MB RAM per task
- **LR Calculation**: ~1.5 hours, ~400 MB RAM per task
- **Combined LR**: ~1 minute, ~95 MB RAM per task
- **Analysis**: ~20 minutes, 96 GB RAM (single task)
- **Prepare Intermediates**: ~10 minutes, 48 GB RAM
- **Plotting**: ~15–20 minutes total, interactive
- **Statistical Tests**: ~30–60 minutes, 64 GB RAM
- **Statistical Report**: ~5 minutes, 10 GB RAM

**Storage Requirements**:
- Simulation files: ~1 GB (1000 files × 1 MB)
- LR files: ~16 GB (1000 files × 16 MB)
- Combined LR files: ~3.3 GB (1000 files × 3.3 MB)
- Analysis outputs: ~7.6 GB
- Total: ~27.3 GB for complete pipeline

**Timeline** (with parallelization):
- Step 1 (Simulation): ~3 hour wall time
- Step 2 (LR Calc): ~1.5 hours wall time
- Step 3 (Combined LR): ~1 minute wall time
- Step 4 (Analysis): ~20 minutes
- Step 4.5 (Prepare Intermediates): ~10 minutes
- Step 5 (Plotting): ~15–20 minutes (interactive)
- Step 6 (Statistical Tests): ~30–60 minutes
- Step 7 (Report): ~10 minutes
- **Total**: ~6 hours from start to final outputs

---

**Last Updated**: April 10, 2026
