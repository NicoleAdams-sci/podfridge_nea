#!/bin/bash
#
# build_pair_ranking_manifest.sh
#
# Scans output/ for pairs_<pop>_<relationship>_n*_chunk*_*.csv files
# restricted to one population and a fixed set of true relationships,
# matches each to its output/combined_LR/combined_LR_...csv counterpart,
# and resolves the single shared unrelated pool file for that population.
# Writes everything to a manifest CSV that submit_pair_ranking.sh reads
# by SLURM_ARRAY_TASK_ID, instead of one file being hardcoded per job.
#
# Usage:
#   code/make_focal_ranking_manifest.sh [POPULATION] [MANIFEST_OUT] [MAX_CHUNK]
#
# MAX_CHUNK (optional) caps which chunk files get included, e.g. MAX_CHUNK=10
# includes chunk1..chunk10 and skips chunk11+. Leave unset/empty to include
# every chunk file found (previous behavior, unchanged).
#
# Example (only need 10k replicates/relationship = 10 chunks of n1000):
#   code/make_focal_ranking_manifest.sh all output/focal_ranking_test/manifest_all_20260702.csv 10

set -euo pipefail

POPULATION=${1:-all}
MANIFEST_OUT=${2:-output/focal_ranking_test/manifest_${POPULATION}_$(date +%Y%m%d_%H%M%S).csv}
MAX_CHUNK=${3:-}

# True-relative types simulated for the focal test (each gets ranked
# against both tested hypotheses inside run_pair_ranking_chunk.R)
RELATIONSHIPS=("parent_child" "full_siblings" "half_siblings")

PAIR_DIR="output"
COMBINED_LR_DIR="output/combined_LR"
UNRELATED_POOL_DIR="output/unrelated_pool"

mkdir -p "$(dirname "$MANIFEST_OUT")"

# -------------------------------------------------------------------------
# Resolve the single unrelated pool file for this population.
# Expected pattern: <pop>_N<size>_combined_unrelated_<timestamp>.csv
# If more than one matches, take the most recently modified.
# -------------------------------------------------------------------------
UNRELATED_POOL_FILE=$(ls -t "${UNRELATED_POOL_DIR}/${POPULATION}"_N*_combined_unrelated_*.csv 2>/dev/null | head -n 1 || true)

if [ -z "$UNRELATED_POOL_FILE" ]; then
  echo "ERROR: no unrelated pool file found matching ${UNRELATED_POOL_DIR}/${POPULATION}_N*_combined_unrelated_*.csv" >&2
  exit 1
fi

echo "Using unrelated pool file: ${UNRELATED_POOL_FILE}"

# -------------------------------------------------------------------------
# Build manifest: one row per pair_file that has a matching combined_LR file
# -------------------------------------------------------------------------
echo "pair_file,combined_lr_file,unrelated_pool_file" > "$MANIFEST_OUT"

n_found=0
n_missing_lr=0
n_skipped_chunk=0

for REL in "${RELATIONSHIPS[@]}"; do
  for PAIR_FILE in "${PAIR_DIR}"/pairs_"${POPULATION}"_"${REL}"_n*_chunk*_*.csv; do
    [ -e "$PAIR_FILE" ] || continue   # glob matched nothing for this relationship

    BASENAME=$(basename "$PAIR_FILE")

    # Extract the chunk number (e.g. "..._chunk10_..." -> 10) for numeric
    # comparison, since chunk10 < chunk2 lexicographically but not numerically.
    CHUNK_NUM=$(echo "$BASENAME" | sed -E 's/.*_chunk([0-9]+)_.*/\1/')

    if [ -n "$MAX_CHUNK" ] && [ "$CHUNK_NUM" -gt "$MAX_CHUNK" ] 2>/dev/null; then
      n_skipped_chunk=$((n_skipped_chunk + 1))
      continue
    fi

    COMBINED_LR_FILE="${COMBINED_LR_DIR}/${BASENAME/pairs_/combined_LR_}"

    if [ -f "$COMBINED_LR_FILE" ]; then
      echo "${PAIR_FILE},${COMBINED_LR_FILE},${UNRELATED_POOL_FILE}" >> "$MANIFEST_OUT"
      n_found=$((n_found + 1))
    else
      echo "WARNING: missing combined_LR file for ${PAIR_FILE} (expected ${COMBINED_LR_FILE})" >&2
      n_missing_lr=$((n_missing_lr + 1))
    fi
  done
done

echo "============================================================="
echo "  Manifest written to : ${MANIFEST_OUT}"
echo "  Population           : ${POPULATION}"
echo "  Max chunk            : ${MAX_CHUNK:-<no limit>}"
echo "  Pair/combined_LR rows: ${n_found}"
echo "  Skipped (missing LR) : ${n_missing_lr}"
echo "  Skipped (> max chunk): ${n_skipped_chunk}"
echo "============================================================="
