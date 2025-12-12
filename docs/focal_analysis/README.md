# Focal Individual Analysis Notes

This directory tracks the lightweight inputs and procedures for focal-individual simulations so we can reproduce focal test runs without reverse engineering ad-hoc scripts.

## Reference configuration

The canonical sandbox config lives at `data/focal_family_configs.csv`. Each row defines a scenario with:

- `scenario_id`: short numeric key for logging and filenames.
- `structure_label`: mnemonic for the family template (used as part of filenames/batch IDs).
- `population`: one of AfAm/Cauc/Hispanic/Asian/all (must exist in `df_allelefreq_combined.csv`).
- `n_focal`: number of focal individuals to generate for that scenario.
- Relationship count columns (`parent_child`, `full_siblings`, `half_siblings`, `cousins`, `second_cousins`, `unrelated`). Leave zero if a relationship is not required.
- `notes`: free text reminder of why that scenario exists.

## Loading and running the scenarios

```r
library(data.table)
library(purrr)
source("code/module7_single_pop_focal_generator.R")
source("code/LR_kinship_utility_functions.R")

allele_frequency_data <- fread("data/df_allelefreq_combined.csv")
allele_frequency_data$frequency[allele_frequency_data$frequency == 0] <- FALLBACK_FREQ
kinship_coefficients <- fread("data/kinship_coefficients.csv")

configs <- fread("data/focal_family_configs.csv")
configs[, structure_label := fifelse(structure_label == "", paste0("scenario_", scenario_id), structure_label)]

split(configs, configs$scenario_id) |>
  purrr::iwalk(function(cfg, id) {
    rel_counts <- as.list(cfg[1, .(parent_child, full_siblings, half_siblings, cousins, second_cousins, unrelated)])
    rel_counts <- rel_counts[unlist(rel_counts) > 0]

    generate_single_pop_focal(
      population = cfg$population[1],
      n_focal = cfg$n_focal[1],
      relationship_counts = rel_counts,
      loci_list = NULL,
      allele_frequency_data = allele_frequency_data,
      kinship_coefficients = kinship_coefficients,
      output_dir = file.path("output", "focal_database", cfg$structure_label[1]),
      custom_datetime = format(Sys.time(), "%Y%m%d_%H%M%S")
    )
  })
```

Notes:
- `relationship_counts` must only include non-zero entries so Module 7 skips unused types cleanly.
- When `loci_list = NULL`, Module 7 simulates every locus present in the allele-frequency table. Pass an explicit list (e.g., `loci_lists$expanded_20`) if you want a restricted panel.
- Set `custom_datetime` when you want to keep outputs synchronized across scenarios.

## Downstream LR testing

1. Stitch the Module 7 output CSVs (they follow the `*_focal_<structure>_<timestamp>.csv` pattern) into the focal database folder under `output/focal_database/<population>/` so the legacy focal LR scripts continue to work.
2. Feed those CSVs into `code/lr_wrapper.R` or the focal-specific LR scripts in `code_graveyard/focal_individual_ranking/` if you are still validating their parity.
3. Combined LRs: reuse `code/combined_lr_wrapper.R` the same way as standard pair simulations—the focal files already expose `focal_A1`, `focal_A2`, etc., so no format bridge is required.

## Tracking experiments

Document every focal run (date, scenario IDs, output paths, LR jobs submitted) below so we can audit assumptions later.

| Date | Scenario IDs | Populations | Output folder | Notes |
|------|--------------|-------------|---------------|-------|
|      |              |             |               |       |
