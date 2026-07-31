# PODFRIDGE Focal Ranking Test

Supplementary workflow for the focal (database-search) ranking test. This
document covers **only** the focal-test-specific scripts. For the core
simulation and likelihood-ratio machinery it builds on — modules 1–5, module 8,
`LR_kinship_utility_functions.R`, and the origin of the `pairs_*.csv` and
`combined_LR_*.csv` inputs — see the main pipeline README,
[`SimulationModules_summary.md`](SimulationModules_summary.md).

---

## What this test asks

Take one member of an already-simulated pair as a "focal" query profile
(individual 1) and its true relative (individual 2). Drop the true relative
into a large database of unrelated individuals, score every database member
against the focal profile by combined LR, rank them, and ask:

> Does the true relative land in the top N (200 / 100 / 50 / 10)?

and how that recovery changes with:

- **true relationship** — parent–child, full siblings, half siblings
- **tested hypothesis** — each pair is scored under both `parent_child` and
  `full_siblings`, regardless of what it truly is
- **loci set** — Core 13 vs Expanded 20
- **database size** — e.g. 10k / 20k / 100k unrelated individuals

Each existing pair is one focal replicate, so 10 chunks × 1,000 pairs =
**10,000 replicates per true relationship**.

---

## Relationship to the main pipeline

The focal test **reuses** work from the main pipeline rather than resimulating:

- **Focal + true relative** come from the existing `pairs_*.csv` files
  (individual 1 = focal, individual 2 = the true relative).
- **The true relative's combined LR** is read straight from the existing
  `output/combined_LR/combined_LR_*.csv` files — not recomputed.
- **Core modules 1–5, module 8, and `LR_kinship_utility_functions.R`** are
  shared, unchanged, and documented in the main README.

The focal layer adds only: database generation, a fast focal-vs-database LR
path, the ranking/outcome recorder, the SLURM orchestration, and the plots.

> **There is no Module 10.** A "Module 10" would have assembled a focal-vs-database
> table. Because this design reuses existing pairs instead of simulating fresh
> focal families, that assembly is done by
> `assemble_unrelated_database_from_existing_pair()` in
> `focal_test_helper_fns.R`.

---

## Files

All scripts live in `code/`. Shared data (`kinship_coefficients.csv`,
`core_CODIS_loci.csv`, `df_allelefreq_combined.csv`) lives in `data/` and is
described in the main README.

### Database generation
| File | Role |
|------|------|
| `generate_unrelated_pool.R` | Build the unrelated database from `all` / `single` / `mixed-proportions` / `mixed-counts` population frequencies (wraps Module 8). Relabels IDs globally so mixed pools don't collide. |
| `generate_unrelated_pool.sh` | SLURM wrapper for the above. |

### Ranking engine (all sourced by `run_focal_ranking.R`)
| File | Role |
|------|------|
| `focal_test_helper_fns.R` | Does the "Module 10" job: reshape the existing true-relative combined LR into the format Module 12 expects, and assemble the focal genotype against every unrelated candidate. |
| `module11_ranking_lr_calculator_fast.R` | Combined LR for every focal-vs-candidate pair, under each tested hypothesis × loci set (wraps the fast single-locus LR + Module 5). |
| `module4_single_locus_LR_fast.R` | Vectorized single-locus LR — one `data.table` join instead of a per-row filter over the allele-frequency table. Drop-in for `module4_single_locus_LR.R` inside Module 11 (~10–100× faster); leaves the original module 4 untouched. |
| `module12_ranking_outcome_recorder.R` | Rank candidates by combined LR within each replicate, find the true relative's rank, record the top 10/50/100/200 flags. |

### Orchestration
| File | Role |
|------|------|
| `make_focal_ranking_manifest.sh` | Match each `pairs_*.csv` to its `combined_LR_*.csv` and an **explicitly chosen** pool file; write a manifest (one row per pair file). Errors rather than guessing if a `--pool-size` matches zero or several pool files. |
| `submit_focal_ranking.sh` | SLURM array driver. Maps each array task to (manifest row, within-file pair chunk) and calls `run_focal_ranking.R`. |
| `run_focal_ranking.R` | Rank one chunk of pair IDs against the pool and write a per-chunk outcomes CSV. |
| `combine_ranking_outcomes.sh` | Concatenate per-chunk `ranking_outcomes_*.csv` (and any `ranking_failures_*.csv`) into one combined CSV per database size. |

### Figures & tables
| File | Role |
|------|------|
| `plot_ranking_summary.R` | **Per-database** figures: recovery bars in top N, rank distribution, optional tied-group size + a summary table. Run once per database size. |
| `compare_ranking_across_databases.R` | **Cross-database** figures: recovery curves (empirical CDF of rank) and per-threshold degradation with Wilson 95% CIs, plus summary tables. Run once over all sizes. |

### Dependency chain
```
run_focal_ranking.R
├── module11_ranking_lr_calculator_fast.R
│   ├── LR_kinship_utility_functions.R        (main pipeline)
│   ├── module1 / module2 / module3           (main pipeline)
│   ├── module4_single_locus_LR.R             (main pipeline, original)
│   ├── module4_single_locus_LR_fast.R        (focal — vectorized swap)
│   └── module5_combined_LR.R                 (main pipeline)
├── module12_ranking_outcome_recorder.R       (focal)
└── focal_test_helper_fns.R                   (focal — the "Module 10" job)
```

---

## Workflow

Steps 1–5 run **once per database size**. Step 6 runs **once**, across all of
them. Substitute your size (e.g. `20000`) throughout.

```bash
# ============================================================================
# PER-SIZE (repeat for each N: 10000, 20000, 100000, ...)
# ============================================================================

# 1. Generate the unrelated pool (once per N)
sbatch code/generate_unrelated_pool.sh all 20000
#   -> output/unrelated_pool/all_N20000_combined_unrelated_<dt>.csv

# 2. Build the manifest (matches pair files + combined_LR files + this pool).
#    Pool is chosen by size; use --pool-file to disambiguate two pools of the
#    same size. --max-chunk 10 keeps chunks 1-10 = 10k replicates/relationship.
code/make_focal_ranking_manifest.sh --pool-size 20000 --max-chunk 10
#   -> output/focal_test_20k/manifest_all_N20000_<dt>.csv

# 3. Rank. Array size = (rows in manifest) x CHUNKS_PER_FILE.
#    3 relationships x 10 chunks = 30 manifest rows; CHUNKS_PER_FILE=100
#    (n1000 pairs / CHUNK_SIZE=10) -> 3000 tasks.
sbatch --array=1-3000 code/submit_focal_ranking.sh \
  output/focal_test_20k/manifest_all_N20000_<dt>.csv \
  10 100 \
  output/focal_test_20k/chunks

# 4. Combine per-chunk outcomes into one CSV (directory-scoped).
code/combine_ranking_outcomes.sh \
  output/focal_test_20k/chunks \
  output/focal_test_20k/combined_all_N20000
#   -> output/focal_test_20k/combined_all_N20000_outcomes.csv

# 5. Per-database figures (subtitle now stamps the database size)
Rscript code/plot_ranking_summary.R \
  output/focal_test_20k/combined_all_N20000_outcomes.csv \
  output/focal_test_20k/figures

# ============================================================================
# ONCE, ACROSS ALL SIZES (needs >= 2 outcome files)
# ============================================================================

# 6. Cross-database comparison. Database size is read from the n_database
#    column, so a mislabeled folder can't mis-color a curve.
Rscript code/compare_ranking_across_databases.R \
  --out output/focal_cross_db \
  output/focal_test_10k/combined_all_N10000_outcomes.csv \
  output/focal_test_20k/combined_all_N20000_outcomes.csv \
  output/focal_test_100k/combined_all_N100000_outcomes.csv
```

### Output layout
```
output/
├── unrelated_pool/
│   └── all_N<size>_combined_unrelated_<dt>.csv
├── focal_test_<size>/                 # one dir per database size
│   ├── manifest_all_N<size>_<dt>.csv
│   ├── chunks/ranking_outcomes_*.csv
│   ├── combined_all_N<size>_outcomes.csv
│   └── figures/
└── focal_cross_db/                    # step 6 output
    ├── recovery_curves_across_databases.{png,pdf}
    ├── recovery_by_threshold_across_databases.{png,pdf}
    ├── recovery_by_threshold_across_databases.csv
    └── rank_summary_across_databases.csv
```

---

## Conventions & gotchas

- **Database size is `n_database - 1`.** The `n_database` column counts every
  ranked candidate, which includes the one inserted true relative. The
  unrelated-pool size is therefore `n_database - 1`.

- **`(batch_id, pair_id)` is unique only *within* a true relationship.** The
  pairs were simulated in a parallel SLURM array and `batch_id` has only
  second resolution, so tasks that started in the same second share a
  `batch_id`; `pair_id` (e.g. `c03_001`) encodes the chunk, not the
  relationship, so it collides across relationships. Any unique-replicate count
  must key on `true_relationship` as well — the plotting scripts already do.

- **Ranking under the wrong hypothesis produces a tied-at-zero block.** Under a
  mismatched hypothesis (especially `parent_child` for a true full/half sib),
  most candidates hit a no-sharing locus and collapse to `combined_LR = 0`.
  `rank(ties.method = "average")` gives that whole block one averaged rank
  ≈ `n_database / 2`. Top-N recovery is the primary, interpretable metric;
  median rank in these cells is an averaging artifact, not a position.

- **Pool selection is never guessed.** `make_focal_ranking_manifest.sh` errors
  if a size matches zero or multiple pool files. Use `--pool-file <path>` to
  pick one explicitly. The chosen pool's `database_label` and `pool_size`
  travel in the manifest, and the size is cross-checked against the pool's
  `_composition_summary_` file.

- **Directory naming.** The manifest script defaults to
  `output/focal_test_<size>/`; keep this consistent across the
  `combine_ranking_outcomes.sh` prefix so step 6's file list is unambiguous.
  The `n_database` column is the source of truth regardless of folder name.

- **The fast single-locus LR was validated for `parent_child` and
  `full_siblings` only.** If you extend the ranking to other tested
  hypotheses, re-validate `module4_single_locus_LR_fast.R` against the original
  `module4_single_locus_LR.R` before trusting the results.
