# STR Simulation and Analysis Pipeline (NEA Version)

This pipeline simulates genotype data for STR loci across different populations and relationship types, processes the data, summarizes the outputs, and generates visualizations comparing matched vs. mismatched allele frequency sources.

---

## ⏱ Order of Execution

1. **`script_NEA.sh`** – SLURM array job submission script
2. **`STR_sims_allPopLR_NEA.R`** – Simulates STR genotypes and calculates LRs
3. **`combine_genotypes_byPopRel_NEA.sh`** – Combines output CSVs by population and relationship
4. **`sims_post_analysis_NEA.R`** – Summarizes LRs and calculates statistics
5. **`plots_known_NEA.R`** – Plots results using population-matched frequencies
6. **`plots_known_vs_tested_NEA.R`** – Plots results comparing matched vs. mismatched frequencies

---

## 📂 File Summaries

### `script_NEA.sh`
- **Purpose:** Launches a SLURM array job to run STR simulations.
- **Inputs:** 
  - `$1` = number of related pairs per group
  - `$2` = number of unrelated pairs
- **SLURM Array Design:**
  - `#SBATCH --array=1-10`: Runs 10 parallel tasks.
  - `--cpus-per-task=8`, `--mem=4g`, `--time=15:00`
- **Output:** Runs `STR_sims_allPopLR_NEA.R` per task, creating subdirectories in `output/`.
- **How to Run:**
  ```bash
  sbatch --array 1-5 code/script_NEA.sh 30 60
  ```

### `STR_sims_allPopLR_NEA.R`
- **Purpose:** Simulates genotypes for STR loci under various relationship types and populations. Computes likelihood ratios (LRs) using allele frequencies from multiple populations.
- **Inputs:**
  - Allele frequencies: `data/df_allelefreq_combined.csv`
  - Core loci list: `data/core_CODIS_loci.csv`
  - Command-line args: number of related and unrelated pairs
- **Outputs (per array task):**
  - `sim_processed_genotypes_task<id>.csv`: Raw genotype + LR results
  - `timing_log_task<id>.csv`: Function timing logs
- **How to Run: (submitted within SLURM job)**
  ```bash
  Rscript code/STR_sims_allPopLR_NEA.R 100 200
  ```

### `combine_genotypes_byPopRel_NEA.sh`
- **Purpose:** Aggregates simulation outputs across array jobs by population and relationship type.
- **Inputs:** Directory names from `output/` (e.g., `output/simulation_*`)
- **Outputs:** 
  - `sim_processed_genotypes_<POP>_<REL>_combined.csv` for each population–relationship pair.
- **Notes:** Also removes temporary per-task population files.
- **How to Run: (Maybe a better way to run this??)**
  ```bash
  bash code/combine_genotypes_byPopRel_NEA.sh simulation_20250508*
  ```

### `sims_post_analysis_NEA.R`
- **Purpose:** Performs post-processing on combined genotype files for one target population.
- **Inputs:** 
  - Directory of combined files (from step 3)
  - Target population (e.g., `AfAm`)
- **Outputs:**
  - `sim_summary_genotypes_<POP>.csv`: Combined LR products by loci set
  - `sim_cutoffs_<POP>.csv`: LR cutoffs at 1%, 0.1%, 0.01% FPR
  - `sim_proportions_exceeding_cutoffs_<POP>.csv`: Proportion of related LRs exceeding cutoffs
  - `sim_lr_summary_stats_<POP>.csv`: LR summary stats
  - `sim_lr_ratio_summary_<POP>.csv`: LR ratios (wrong/correct)
  - `timing_log_summary_<POP>.csv`: Processing times
- **How to Run:**
  ```bash
  # For one population
  Rscript code/sims_post_analysis_NEA.R output AfAm output/AfAm_summary
  
  # For all populations
  for POP in AfAm Cauc Asian Hispanic; do echo $POP; 
  Rscript code/sims_post_analysis_NEA.R output $POP output/$POP\_summary;
  done
  ```

### `plots_known_NEA.R`
- **Purpose:** Generates plots **only using matched allele frequency sources**.
- **Input:** Folder with summary outputs (from previous step)
- **Output Plots:**
  - Boxplots and violin plots of LRs by relationship and loci set
  - Mean LR line plots
  - Proportion bar charts and heatmaps (cutoff exceedance)
- **How to Run:**
  ```bash
  Rscript code/plots_known_NEA.R output test_plots_known
  ```

### `plots_known_vs_tested_NEA.R`
- **Purpose:** Generates plots comparing **matched and mismatched allele frequency** effects.
- **Input:** Same summary folders
- **Output Plots:**
  - Boxplots and violin plots by mismatched frequency source
  - Ratio heatmaps and boxplots (wrong LR / correct LR)
- **How to Run:**
  ```bash
  Rscript code/plots_known_vs_tested_NEA.R output test_plots_knownTested
  ```