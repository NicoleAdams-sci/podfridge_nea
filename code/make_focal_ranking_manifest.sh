#!/bin/bash
#
# make_focal_ranking_manifest.sh
#
# Scans output/ for pairs_<pop>_<relationship>_n*_chunk*_*.csv files
# restricted to one population and a fixed set of true relationships,
# matches each to its output/combined_LR/combined_LR_...csv counterpart,
# and pairs them with an EXPLICITLY CHOSEN unrelated pool file. Writes
# everything to a manifest CSV that submit_focal_ranking.sh reads by
# SLURM_ARRAY_TASK_ID.
#
# The unrelated pool is selected by size (--pool-size) or by exact path
# (--pool-file). It is never guessed from modification time: if a size
# matches zero or multiple files, this script errors out rather than
# silently picking one.
#
# Usage:
#   code/make_focal_ranking_manifest.sh --pool-size 20000 [options]
#   code/make_focal_ranking_manifest.sh --pool-file <path> [options]
#
# Options:
#   --pool-size N     Select output/unrelated_pool/<pop>_N<N>_combined_unrelated_*.csv
#   --pool-file PATH  Use this exact pool file (overrides --pool-size)
#   --population POP  Population prefix for pair files and pool. Default: all
#   --max-chunk N     Only include chunk1..chunkN of the pair files. Default: no limit
#   --out PATH        Manifest output path.
#                     Default: output/focal_test_<size>/manifest_<pop>_N<size>_<datetime>.csv
#   --list            List the available pool files and exit
#   -h, --help        Show this message
#
# Examples:
#   # 20k database, 10 chunks of n1000 pairs = 10k replicates per relationship
#   code/make_focal_ranking_manifest.sh --pool-size 20000 --max-chunk 10
#
#   # 100k database, explicit output path
#   code/make_focal_ranking_manifest.sh --pool-size 100000 --max-chunk 10 \
#     --out output/focal_test_100k/manifest_all_N100000.csv
#
#   # Disambiguate two pools of the same size
#   code/make_focal_ranking_manifest.sh \
#     --pool-file output/unrelated_pool/all_N20000_combined_unrelated_20260716_120513.csv
#
# Manifest columns:
#   pair_file,combined_lr_file,unrelated_pool_file,database_label,pool_size
# The first three are what submit_focal_ranking.sh cuts out by position;
# the last two exist so the database identity travels with the manifest
# instead of living only in a directory name.

set -euo pipefail

PAIR_DIR="output"
COMBINED_LR_DIR="output/combined_LR"
UNRELATED_POOL_DIR="output/unrelated_pool"

# True-relative types simulated for the focal test (each gets ranked
# against both tested hypotheses inside run_focal_ranking.R)
RELATIONSHIPS=("parent_child" "full_siblings" "half_siblings")

usage() {
  sed -n '2,50p' "$0" | sed 's/^# \?//'
}

list_pools() {
  echo "Available combined unrelated pool files in ${UNRELATED_POOL_DIR}:"
  echo
  local found=0
  shopt -s nullglob
  for f in "${UNRELATED_POOL_DIR}"/*_combined_unrelated_*.csv; do
    printf '  %-70s  %s\n' "$(basename "$f")" "$(date -r "$f" '+%Y-%m-%d %H:%M')"
    found=1
  done
  shopt -u nullglob
  [ "$found" -eq 0 ] && echo "  (none found)"
  echo
}

# -------------------------------------------------------------------------
# Parse arguments
# -------------------------------------------------------------------------

POPULATION="all"
POOL_SIZE=""
POOL_FILE=""
MANIFEST_OUT=""
MAX_CHUNK=""

if [ "$#" -eq 0 ]; then
  usage >&2
  echo "ERROR: one of --pool-size or --pool-file is required." >&2
  exit 1
fi

while [ "$#" -gt 0 ]; do
  case "$1" in
    --pool-size)   POOL_SIZE="${2:?--pool-size needs a value}";   shift 2 ;;
    --pool-file)   POOL_FILE="${2:?--pool-file needs a value}";   shift 2 ;;
    --population)  POPULATION="${2:?--population needs a value}"; shift 2 ;;
    --max-chunk)   MAX_CHUNK="${2:?--max-chunk needs a value}";   shift 2 ;;
    --out)         MANIFEST_OUT="${2:?--out needs a value}";      shift 2 ;;
    --list)        list_pools; exit 0 ;;
    -h|--help)     usage; exit 0 ;;
    *)
      echo "ERROR: unknown argument: $1" >&2
      echo "Run with --help for usage." >&2
      exit 1
      ;;
  esac
done

if [ -n "$POOL_SIZE" ] && [ -n "$POOL_FILE" ]; then
  echo "ERROR: --pool-size and --pool-file are mutually exclusive. Pick one." >&2
  exit 1
fi

if [ -z "$POOL_SIZE" ] && [ -z "$POOL_FILE" ]; then
  echo "ERROR: one of --pool-size or --pool-file is required." >&2
  echo >&2
  list_pools >&2
  exit 1
fi

if [ -n "$POOL_SIZE" ] && ! [[ "$POOL_SIZE" =~ ^[0-9]+$ ]]; then
  echo "ERROR: --pool-size must be a positive integer (got '${POOL_SIZE}')." >&2
  exit 1
fi

if [ -n "$MAX_CHUNK" ] && ! [[ "$MAX_CHUNK" =~ ^[0-9]+$ ]]; then
  echo "ERROR: --max-chunk must be a positive integer (got '${MAX_CHUNK}')." >&2
  exit 1
fi

# -------------------------------------------------------------------------
# Resolve the unrelated pool file
#
# The glob is anchored on _N<size>_combined_unrelated_ rather than _N<size>*
# because all_N1000* would otherwise also match all_N10000 and all_N100000.
# -------------------------------------------------------------------------

if [ -n "$POOL_FILE" ]; then

  if [ ! -f "$POOL_FILE" ]; then
    echo "ERROR: --pool-file does not exist: ${POOL_FILE}" >&2
    exit 1
  fi

else

  shopt -s nullglob
  matches=("${UNRELATED_POOL_DIR}/${POPULATION}_N${POOL_SIZE}_combined_unrelated_"*.csv)
  shopt -u nullglob

  if [ "${#matches[@]}" -eq 0 ]; then
    echo "ERROR: no pool file matching ${UNRELATED_POOL_DIR}/${POPULATION}_N${POOL_SIZE}_combined_unrelated_*.csv" >&2
    echo >&2
    echo "Generate it first:" >&2
    echo "  sbatch code/generate_unrelated_pool.sh all ${POOL_SIZE}" >&2
    echo >&2
    list_pools >&2
    exit 1
  fi

  if [ "${#matches[@]}" -gt 1 ]; then
    echo "ERROR: ${#matches[@]} pool files match population=${POPULATION}, size=${POOL_SIZE}:" >&2
    printf '  %s\n' "${matches[@]}" >&2
    echo >&2
    echo "Refusing to guess. Pass the one you want explicitly with --pool-file." >&2
    exit 1
  fi

  POOL_FILE="${matches[0]}"

fi

# Derive database_label ("all_N20000") by stripping the _combined_unrelated_<datetime>.csv suffix
POOL_BASENAME=$(basename "$POOL_FILE")
DATABASE_LABEL=$(echo "$POOL_BASENAME" | sed -E 's/_combined_unrelated_[0-9]{8}_[0-9]{6}\.csv$//')

# Recover the size from the label when it came in via --pool-file
if [ -z "$POOL_SIZE" ]; then
  POOL_SIZE=$(echo "$DATABASE_LABEL" | sed -nE 's/.*_N([0-9]+)$/\1/p')
  if [ -z "$POOL_SIZE" ]; then
    echo "WARNING: could not parse a size out of '${DATABASE_LABEL}'." >&2
    echo "         pool_size will be recorded as NA in the manifest." >&2
    POOL_SIZE="NA"
  fi
fi

# Cross-check against the composition summary written alongside the pool,
# which records the size that was actually requested.
SUMMARY_FILE="${POOL_FILE/_combined_unrelated_/_composition_summary_}"
POOL_SIZE_CONFIRMED="(no composition summary found)"
if [ -f "$SUMMARY_FILE" ] && [ "$POOL_SIZE" != "NA" ]; then
  SUMMARY_N=$(awk -F',' 'NR > 1 { total += $4 } END { print total+0 }' "$SUMMARY_FILE")
  if [ "$SUMMARY_N" = "$POOL_SIZE" ]; then
    POOL_SIZE_CONFIRMED="confirmed (${SUMMARY_N} individuals)"
  else
    echo "WARNING: filename says N=${POOL_SIZE} but ${SUMMARY_FILE}" >&2
    echo "         sums to ${SUMMARY_N} individuals. Check the pool before running." >&2
    POOL_SIZE_CONFIRMED="MISMATCH (summary says ${SUMMARY_N})"
  fi
fi

# -------------------------------------------------------------------------
# Default manifest path, e.g. output/focal_test_20k/manifest_all_N20000_<dt>.csv
# -------------------------------------------------------------------------

if [ -z "$MANIFEST_OUT" ]; then
  if [ "$POOL_SIZE" != "NA" ] && [ $((POOL_SIZE % 1000)) -eq 0 ] && [ "$POOL_SIZE" -ge 1000 ]; then
    SIZE_TAG="$((POOL_SIZE / 1000))k"
  else
    SIZE_TAG="N${POOL_SIZE}"
  fi
  MANIFEST_OUT="output/focal_test_${SIZE_TAG}/manifest_${POPULATION}_N${POOL_SIZE}_$(date +%Y%m%d_%H%M%S).csv"
fi

mkdir -p "$(dirname "$MANIFEST_OUT")"

echo "Using unrelated pool file : ${POOL_FILE}"
echo "Database label            : ${DATABASE_LABEL}"
echo "Pool size                 : ${POOL_SIZE} - ${POOL_SIZE_CONFIRMED}"
echo

# -------------------------------------------------------------------------
# Build manifest: one row per pair_file that has a matching combined_LR file
# -------------------------------------------------------------------------

echo "pair_file,combined_lr_file,unrelated_pool_file,database_label,pool_size" > "$MANIFEST_OUT"

n_found=0
n_missing_lr=0
n_skipped_chunk=0

for REL in "${RELATIONSHIPS[@]}"; do
  shopt -s nullglob
  for PAIR_FILE in "${PAIR_DIR}"/pairs_"${POPULATION}"_"${REL}"_n*_chunk*_*.csv; do

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
      echo "${PAIR_FILE},${COMBINED_LR_FILE},${POOL_FILE},${DATABASE_LABEL},${POOL_SIZE}" >> "$MANIFEST_OUT"
      n_found=$((n_found + 1))
    else
      echo "WARNING: missing combined_LR file for ${PAIR_FILE} (expected ${COMBINED_LR_FILE})" >&2
      n_missing_lr=$((n_missing_lr + 1))
    fi
  done
  shopt -u nullglob
done

if [ "$n_found" -eq 0 ]; then
  echo "ERROR: manifest is empty - no pair/combined_LR file pairs found for population '${POPULATION}'." >&2
  exit 1
fi

echo "============================================================="
echo "  Manifest written to  : ${MANIFEST_OUT}"
echo "  Population           : ${POPULATION}"
echo "  Database label       : ${DATABASE_LABEL}"
echo "  Pool file            : ${POOL_FILE}"
echo "  Max chunk            : ${MAX_CHUNK:-<no limit>}"
echo "  Pair/combined_LR rows: ${n_found}"
echo "  Skipped (missing LR) : ${n_missing_lr}"
echo "  Skipped (> max chunk): ${n_skipped_chunk}"
echo "============================================================="
echo
echo "Next step - submit the ranking array. With CHUNKS_PER_FILE=100:"
echo
echo "  sbatch --array=1-$((n_found * 100)) code/submit_focal_ranking.sh \\"
echo "    ${MANIFEST_OUT} \\"
echo "    10 100 \\"
echo "    $(dirname "$MANIFEST_OUT")/chunks"
echo
