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
# Example:
#   code/combine_ranking_outcomes.sh \
#     output/focal_ranking_test/ \
#     output/focal_ranking_test/combined_all_20260709
#
# Produces:
#   <PREFIX>_outcomes.csv
#   <PREFIX>_failures.csv   (only if any ranking_failures_*.csv files exist)

set -euo pipefail

OUTPUT_DIR=${1:?"Usage: code/combine_ranking_outcomes.sh <OUTPUT_DIR> [COMBINED_OUT_PREFIX]"}
PREFIX=${2:-${OUTPUT_DIR%/}/combined_ranking_$(date +%Y%m%d_%H%M%S)}

OUTCOMES_OUT="${PREFIX}_outcomes.csv"
FAILURES_OUT="${PREFIX}_failures.csv"

shopt -s nullglob
outcome_files=("${OUTPUT_DIR%/}"/ranking_outcomes_*.csv)
failure_files=("${OUTPUT_DIR%/}"/ranking_failures_*.csv)

if [ ${#outcome_files[@]} -eq 0 ]; then
  echo "ERROR: no ranking_outcomes_*.csv files found in ${OUTPUT_DIR}" >&2
  exit 1
fi

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
