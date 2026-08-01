#!/bin/bash
#
# combine_ranking_outcomes.sh
#
# Concatenates all per-chunk ranking_outcomes_*.csv (and, if present,
# ranking_failures_*.csv) files from a focal ranking array run into single
# combined CSVs, keeping one header row instead of one per chunk file.
#
# Usage:
#   code/combine_ranking_outcomes.sh <OUTPUT_DIR> [COMBINED_OUT_PREFIX]
#
# OUTPUT_DIR is the directory holding the per-chunk ranking_outcomes_*.csv
# (normally output/focal_test_<size>/chunks).
#
# If COMBINED_OUT_PREFIX is omitted, the output is auto-named
#   <parent-of-OUTPUT_DIR>/combined_all_N<poolsize>_outcomes.csv
# where <poolsize> is read from the data (the n_database column of the first
# outcome file, minus the 1 inserted true relative). This keeps every
# database's combined file on the same convention as the pool files
# (all_N<size>_combined_unrelated_...) and manifests (manifest_all_N<size>_...).
#
# Pass COMBINED_OUT_PREFIX explicitly to override the auto-name.
#
# Examples:
#   # Auto-named -> output/focal_test_20k/combined_all_N20000_outcomes.csv
#   code/combine_ranking_outcomes.sh output/focal_test_20k/chunks
#
#   # Explicit prefix
#   code/combine_ranking_outcomes.sh output/focal_test_20k/chunks \
#     output/focal_test_20k/combined_all_N20000
#
# Produces:
#   <PREFIX>_outcomes.csv
#   <PREFIX>_failures.csv   (only if any ranking_failures_*.csv files exist)

set -euo pipefail

OUTPUT_DIR=${1:?"Usage: code/combine_ranking_outcomes.sh <OUTPUT_DIR> [COMBINED_OUT_PREFIX]"}
EXPLICIT_PREFIX=${2:-}

shopt -s nullglob
outcome_files=("${OUTPUT_DIR%/}"/ranking_outcomes_*.csv)
failure_files=("${OUTPUT_DIR%/}"/ranking_failures_*.csv)

if [ ${#outcome_files[@]} -eq 0 ]; then
  echo "ERROR: no ranking_outcomes_*.csv files found in ${OUTPUT_DIR}" >&2
  exit 1
fi

# -------------------------------------------------------------------------
# Determine the output prefix.
# Explicit arg wins; otherwise auto-name from the data so every database's
# combined file follows the same convention: combined_all_N<poolsize>.
# -------------------------------------------------------------------------
if [ -n "$EXPLICIT_PREFIX" ]; then
  PREFIX="$EXPLICIT_PREFIX"
else
  first_file="${outcome_files[0]}"

  # Read n_database from the first data row, locating the column by header
  # name so this is robust to added columns (e.g. tied_group_size).
  N_DB=$(awk -F',' '
    NR==1 { for (i = 1; i <= NF; i++) if ($i == "n_database") c = i; next }
    NR==2 { if (c) print $c; exit }
  ' "$first_file")

  PARENT_DIR=$(dirname "${OUTPUT_DIR%/}")

  if [[ "$N_DB" =~ ^[0-9]+$ ]] && [ "$N_DB" -gt 1 ]; then
    POOL_SIZE=$((N_DB - 1))
    PREFIX="${PARENT_DIR}/combined_all_N${POOL_SIZE}"
    echo "Auto-named from n_database=${N_DB} (pool size ${POOL_SIZE}) -> ${PREFIX}_outcomes.csv"
  else
    PREFIX="${PARENT_DIR}/combined_all_$(date +%Y%m%d_%H%M%S)"
    echo "WARNING: could not read n_database from ${first_file}; using dated name ${PREFIX}." >&2
  fi
fi

OUTCOMES_OUT="${PREFIX}_outcomes.csv"
FAILURES_OUT="${PREFIX}_failures.csv"

# Header from the first file, then data rows (skip each file's own header) from all
head -n 1 "${outcome_files[0]}" > "$OUTCOMES_OUT"
for f in "${outcome_files[@]}"; do
  tail -n +2 "$f" >> "$OUTCOMES_OUT"
done

n_rows=$(($(wc -l < "$OUTCOMES_OUT") - 1))
echo "Combined ${#outcome_files[@]} outcome files -> ${OUTCOMES_OUT} (${n_rows} data rows)"

if [ ${#failure_files[@]} -gt 0 ]; then
  head -n 1 "${failure_files[0]}" > "$FAILURES_OUT"
  for f in "${failure_files[@]}"; do
    tail -n +2 "$f" >> "$FAILURES_OUT"
  done
  n_fail=$(($(wc -l < "$FAILURES_OUT") - 1))
  echo "Combined ${#failure_files[@]} failure files -> ${FAILURES_OUT} (${n_fail} failed pairs) - worth checking these"
else
  echo "No ranking_failures_*.csv files found - no failed pairs to combine."
fi
