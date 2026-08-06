# podfridge_nea

STR (Short Tandem Repeat) profile simulation and kinship likelihood ratio calculation pipeline for forensic genetics research.

## Overview

This repository contains a modular R pipeline for simulating STR genetic profiles and calculating likelihood ratios (LRs) for kinship testing. The pipeline enables evaluation of how population-specific allele frequencies affect kinship hypothesis testing across different relationship types.

### Where to look

This file is the map. The detailed documentation lives in two places:

| Document | Covers |
|---|---|
| [`SimulationModules_summary.md`](SimulationModules_summary.md) | **The main reference.** Every module and wrapper in depth, the full start-to-finish workflow with SLURM commands and resource requirements, output file schemas, and the test suite. |
| [`README_focal_ranking_test.md`](README_focal_ranking_test.md) | The supplementary focal (database-search) ranking test: does a true relative rank in the top N of a database of unrelated individuals? |

## Repository Structure

```
podfridge_nea/
├── code/                  # All R modules, scripts, and wrappers
├── data/                  # Input data (allele frequencies, kinship coefficients)
├── analysis/              # R Markdown reports and their rendered HTML
├── tests/                 # testthat unit and integration tests
├── docs/                  # Reference tables (e.g. SLURM array assignments)
├── output/                # All simulation and analysis outputs
│   ├── LR/                # Single-locus LR results
│   └── combined_LR/       # Multi-locus combined LR results
├── logs/                  # SLURM job logs
└── code_graveyard/        # Archived/deprecated code
```

---

## Data Files (`data/`)

| File | Description |
|------|-------------|
| `df_allelefreq_combined.csv` | Population-specific allele frequencies for STR loci |
| `kinship_coefficients.csv` | IBD sharing coefficients (k0, k1, k2) for relationship types |
| `core_CODIS_loci.csv` | List of core CODIS STR loci |

---

## Core Modules (`code/`)

The pipeline is built on 9 interconnected R modules:

### Foundation Module

| File | Purpose |
|------|---------|
| `LR_kinship_utility_functions.R` | **Core foundation** — Loads allele frequencies, defines loci sets (core_13, identifiler_15, expanded_20, supplementary, autosomal_29), kinship coefficients, and provides core LR calculation functions (`calculate_likelihood_ratio()`, `kinship_calculation()`) |

### Simulation Modules (Modules 1-3)

| File | Purpose |
|------|---------|
| `module1_allele_simulator.R` | Simulates a single allele based on population-specific frequencies via `simulate_allele()` |
| `module2_STR_profile_simulator.R` | Generates complete diploid STR profiles via `simulate_str_profile()` (calls Module 1 twice per locus) |
| `module3_related_individual_simulator.R` | Simulates genetically related individuals based on IBD sharing patterns (k0/k1/k2) via `simulate_related_individual()`, `simulate_individual_pair()`, and `simulate_multiple_pairs()` |

### Likelihood Ratio Calculation Modules (Modules 4-5)

| File | Purpose |
|------|---------|
| `module4_single_locus_LR.R` | Calculates per-locus likelihood ratios for pairs via `calculate_single_locus_lr()` |
| `module5_combined_LR.R` | Combines single-locus LRs into multi-locus LRs by multiplication via `calculate_combined_lr()` for different loci sets |

### Batch Generation Modules (Modules 6-8)

| File | Purpose | Status |
|------|---------|--------|
| `module6_single_combo_pair_generator.R` | Batch generates pairs for specific population-relationship combinations | Unused — superseded by direct `sim_pairs.R` calls |
| `module7_single_pop_focal_generator.R` | Generates focal individuals with flexible family structures | Unused — the focal ranking test reuses existing pairs instead of simulating families. Retained and unit-tested |
| `module8_unrelated_pool_generator.R` | Generates pools of unrelated individuals for null distributions | **In production** — drives `generate_unrelated_pool.R`, step 1 of the focal ranking test |

### Analysis Module (Module 9)

| File | Purpose |
|------|---------|
| `module9_combinedLR_stats_functions.R` | Statistical analysis functions: `calculate_summary_stats()`, `calculate_ratio_stats()`, `calculate_cutoffs()`, `calculate_proportions_exceeding_cutoffs()` |

### Focal Ranking Modules (Modules 10-11)

Supplementary — see [`README_focal_ranking_test.md`](README_focal_ranking_test.md).

| File | Purpose |
|------|---------|
| `focal_test_helper_fns.R` | Assembles the focal profile against every unrelated database candidate, reusing existing pairs |
| `module4_single_locus_LR_fast.R` | Vectorized `data.table`-join rewrite of Module 4 (~10-100x faster). Validated for `parent_child` and `full_siblings` only |
| `module10_ranking_lr_calculator_fast.R` | Combined LR for every focal-vs-candidate pair, per tested hypothesis x loci set |
| `module11_ranking_outcome_recorder.R` | Ranks candidates, finds the true relative's rank, records top 10/50/100/200 flags |

---

## Wrapper Scripts (`code/`)

These scripts orchestrate the modules for batch processing:

| File | Purpose | Usage |
|------|---------|-------|
| `sim_pairs.R` | Simulates pairs of individuals for a given population and relationship | `Rscript code/sim_pairs.R <POP> <RELATIONSHIP> <N_PAIRS> [CHUNK_NUM]` |
| `lr_wrapper.R` | Calculates single-locus LRs for all tested relationships and populations | `Rscript code/lr_wrapper.R <PAIRS_CSV_FILE>` |
| `combined_lr_wrapper.R` | Calculates combined multi-locus LRs across loci sets | `Rscript code/combined_lr_wrapper.R <LR_CSV_FILE>` |
| `analyze_lr_outputs.R` | Aggregates chunked outputs and generates summary statistics | `Rscript code/analyze_lr_outputs.R [OUTPUT_DIR]` |

---

## SLURM Job Scripts (`code/`)

For HPC cluster submission:

| File | Purpose |
|------|---------|
| `sim_pairs.sh` | SLURM array job for parallel pair simulation (5 pops × 6 relationships × 10 chunks) |
| `sim_pairs_unrelated.sh` | SLURM wrapper for simulating unrelated pairs |
| `lr_submission.sh` | Generates file list and submits LR calculation jobs |
| `lr_wrapper.sh` | SLURM wrapper for `lr_wrapper.R` |
| `combined_lr_submission.sh` | Submits combined LR calculation jobs |
| `combined_lr.sh` | SLURM wrapper for `combined_lr_wrapper.R` |
| `analyze_lr_outputs.sh` | SLURM wrapper for `analyze_lr_outputs.R` |
| `prepare_combined_lr_intermediates.sh` | SLURM wrapper for the plotting-intermediates step |
| `analyze_locus_inflation.sh` | SLURM wrapper for per-locus inflation analysis |
| `run_statistical_tests.sh` | SLURM wrapper for the inferential test suite |
| `simulation_analysis.sh` | Workflow script for running analysis pipeline |
| `generate_unrelated_pool.sh` | Focal ranking test — builds an unrelated database |
| `make_focal_ranking_manifest.sh` | Focal ranking test — builds the run manifest |
| `submit_focal_ranking.sh` | Focal ranking test — SLURM array ranking driver |
| `combine_ranking_outcomes.sh` | Focal ranking test — concatenates per-chunk outcomes |

**Note:** The analysis and simulation steps have `.sh` SLURM submission scripts. The publication plotting scripts do **not** — they are run interactively with `Rscript`.

---

## Plotting Scripts (`code/`)

All publication plotting scripts run interactively (`Rscript`), not via SLURM, and share the Okabe-Ito colorblind-safe palette.

| File | Purpose |
|------|---------|
| `plots_matched_publication.R` | LR distributions when population is correctly matched (violin + per-population boxplots) |
| `plots_mismatched_population.R` | Robustness to wrong population frequencies (line plots + inflation heatmaps) |
| `plots_mismatched_relationship.R` | Discrimination between relationship hypotheses |
| `plots_cutoffs_publication.R` | Classification performance at fixed FPR thresholds |
| `plots_locus_inflation.R` | Per-locus LR inflation and heterozygosity |
| `plot_ranking_summary.R` | Focal ranking test — per-database figures |
| `compare_ranking_across_databases.R` | Focal ranking test — cross-database comparison |

> The earlier `plots_matched.R`, `plots_mismatched.R`, and `plots_proportion_exceeding_cutoffs.R` have been superseded; `plots_matched.R` is archived in `code_graveyard/`.

---

## Pipeline Flow

```
┌─────────────────────────────────────────────────────────────────┐
│         LR_kinship_utility_functions.R (Foundation)             │
│    (Allele frequencies, kinship coefficients, core LR logic)    │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│  Module 1 (Allele) → Module 2 (Profile) → Module 3 (Related)    │
│                    SIMULATION LAYER                              │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                     sim_pairs.R / sim_pairs.sh                  │
│              Generate pairs CSV files → output/                  │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│            lr_wrapper.R + Module 4 (Single Locus LR)            │
│           Calculate per-locus LRs → output/LR/                  │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│          combined_lr_wrapper.R + Module 5 (Combined LR)         │
│         Multiply across loci → output/combined_LR/              │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│        analyze_lr_outputs.R + Module 9 (Statistics)             │
│      Aggregate & summarize → output/lr_analysis_<timestamp>/    │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│            plots_*.R (Visualization Scripts)                    │
│                Generate figures for analysis                     │
└─────────────────────────────────────────────────────────────────┘
```

---

## Output Structure

| Directory | Contents |
|-----------|----------|
| `output/` | Raw simulated pairs CSV files (`pairs_<POP>_<REL>_n<N>_chunk<C>_<DATE>.csv`) |
| `output/LR/` | Single-locus LR results (`LR_<POP>_<REL>_n<N>_chunk<C>_<DATE>.csv`) |
| `output/combined_LR/` | Combined multi-locus LR results (`combined_LR_<POP>_<REL>_*.csv`) |
| `output/lr_analysis_<timestamp>/` | Aggregated analysis outputs including: |
| | • `combined_LR_all.rds` — All combined LRs |
| | • `combined_LR_match.csv` — Population-matched results |
| | • `combined_LR_mismatch.csv` — Population-mismatched results |
| | • `combined_LR_summary_stats.csv` — Summary statistics |
| | • `combined_LR_ratio_summary.csv` — LR ratio analysis |

---

## Supported Parameters

### Populations
- `AfAm` (African American)
- `Cauc` (Caucasian)
- `Hispanic`
- `Asian`
- `all` (pooled frequencies)

### Relationships
| Relationship | k0 | k1 | k2 |
|-------------|-----|-----|-----|
| `parent_child` | 0 | 1 | 0 |
| `full_siblings` | 0.25 | 0.5 | 0.25 |
| `half_siblings` | 0.5 | 0.5 | 0 |
| `cousins` | 0.875 | 0.125 | 0 |
| `second_cousins` | 0.9375 | 0.0625 | 0 |
| `unrelated` | 1 | 0 | 0 |

### Loci Sets
- `core_13` — Core 13 CODIS loci
- `identifiler_15` — Identifiler 15 loci
- `expanded_20` — Expanded 20 CODIS loci
- `supplementary` — Additional supplementary loci
- `autosomal_29` — Full autosomal 29 loci panel

---

## Quick Start

### 1. Simulate pairs
```bash
# Simulate 1000 pairs of African American parent-child relationships
Rscript code/sim_pairs.R AfAm parent_child 1000
```

### 2. Calculate single-locus LRs
```bash
Rscript code/lr_wrapper.R output/pairs_AfAm_parent_child_n1000_20251211.csv
```

### 3. Calculate combined LRs
```bash
Rscript code/combined_lr_wrapper.R output/LR/LR_AfAm_parent_child_n1000_20251211.csv
```

### 4. Run analysis and generate plots
```bash
Rscript code/analyze_lr_outputs.R output/my_analysis
Rscript code/prepare_combined_lr_intermediates.R output/my_analysis
Rscript code/plots_matched_publication.R output/my_analysis
```

### 5. Run the test suite
```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```

### HPC Cluster (SLURM)
```bash
# Submit all pair simulations as array job
sbatch code/sim_pairs.sh

# Submit LR calculations
bash code/lr_submission.sh

# Submit combined LR calculations
bash code/combined_lr_submission.sh
```

---

## Dependencies

- R (≥4.0)
- R packages: `dplyr`, `data.table`, `tidyverse`, `tictoc`, `scales`
