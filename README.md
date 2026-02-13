# podfridge_nea

STR (Short Tandem Repeat) profile simulation and kinship likelihood ratio calculation pipeline for forensic genetics research.

## Overview

This repository contains a modular R pipeline for simulating STR genetic profiles and calculating likelihood ratios (LRs) for kinship testing. The pipeline enables evaluation of how population-specific allele frequencies affect kinship hypothesis testing across different relationship types.

## Repository Structure

```
podfridge_nea/
├── code/                  # All R modules, scripts, and wrappers
├── data/                  # Input data (allele frequencies, kinship coefficients)
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

| File | Purpose |
|------|---------|
| `module6_single_combo_pair_generator.R` | Batch generates pairs for specific population-relationship combinations |
| `module7_single_pop_focal_generator.R` | Generates focal individuals with flexible family structures |
| `module8_unrelated_pool_generator.R` | Generates pools of unrelated individuals for null distributions |

### Analysis Module (Module 9)

| File | Purpose |
|------|---------|
| `module9_combinedLR_stats_functions.R` | Statistical analysis functions: `calculate_summary_stats()`, `calculate_ratio_stats()`, `calculate_cutoffs()`, `calculate_proportions_exceeding_cutoffs()` |

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
| `simulation_analysis.sh` | Workflow script for running analysis pipeline |
| `plots_mismatched.sh` | SLURM wrapper for `plots_mismatched.R` |
| `plots_proportion_exceeding_cutoffs.sh` | SLURM wrapper for `plots_proportion_exceeding_cutoffs.R` |

**Note:** Most R wrapper and plotting scripts have corresponding `.sh` SLURM submission scripts for cluster execution.

---

## Plotting Scripts (`code/`)

| File | Purpose |
|------|---------|
| `plots_matched.R` | Plots LR distributions when population is correctly matched |
| `plots_mismatched.R` | Plots LR distributions when population is mismatched |
| `plots_proportion_exceeding_cutoffs.R` | Plots proportion of pairs exceeding LR cutoffs |

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
Rscript code/plots_matched.R output/my_analysis
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
