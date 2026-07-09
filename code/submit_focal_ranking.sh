#!/bin/bash
#SBATCH --job-name=focal_rank
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=logs/focal_rank_%A_%a.out
#SBATCH --error=logs/focal_rank_%A_%a.err

# Usage:
#   sbatch --array=1-$((N_FILES * CHUNKS_PER_FILE)) code/submit_pair_ranking.sh \
#     <MANIFEST_FILE> <CHUNK_SIZE> <CHUNKS_PER_FILE> <OUTPUT_DIR>
#
# MANIFEST_FILE is produced by code/build_pair_ranking_manifest.sh and has
# columns: pair_file,combined_lr_file,unrelated_pool_file (header on line 1,
# one data row per pair_file).
#
# CHUNKS_PER_FILE = ceil(pairs_per_file / CHUNK_SIZE), sized for the
# largest file. e.g. n1000-pair files with CHUNK_SIZE=10 -> 100. Tasks
# that land past the real number of pairs in a smaller file will exit
# with a harmless "beyond available pair IDs" error from the R script.
#
# Example (3 files in manifest, n1000 pairs, chunk_size=10 -> 100/file):
#   sbatch --array=1-300 code/submit_pair_ranking.sh \
#     output/focal_ranking_test/manifest_all_20260702.csv \
#     10 100 \
#     output/focal_ranking_test/pair_rank_chunks

mkdir -p logs
module load Rtidyverse

MANIFEST_FILE=$1
CHUNK_SIZE=$2
CHUNKS_PER_FILE=$3
OUTPUT_DIR=$4

if [ -z "$MANIFEST_FILE" ] || [ -z "$CHUNK_SIZE" ] || [ -z "$CHUNKS_PER_FILE" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Usage: sbatch --array=1-N code/submit_pair_ranking.sh <MANIFEST_FILE> <CHUNK_SIZE> <CHUNKS_PER_FILE> <OUTPUT_DIR>"
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"

# -------------------------------------------------------------------------
# Map this array task to (manifest row, within-file pair-id chunk)
# Same div/mod pattern as sim_pairs.sh's population/relationship mapping.
# -------------------------------------------------------------------------
FILE_INDEX=$(( (SLURM_ARRAY_TASK_ID - 1) / CHUNKS_PER_FILE + 1 ))     # 1-based data row in manifest
CHUNK_INDEX=$(( (SLURM_ARRAY_TASK_ID - 1) % CHUNKS_PER_FILE + 1 ))    # 1-based sub-chunk within that file

MANIFEST_ROW=$(sed -n "$((FILE_INDEX + 1))p" "$MANIFEST_FILE")   # +1 to skip header

if [ -z "$MANIFEST_ROW" ]; then
    echo "No manifest row at index ${FILE_INDEX} (array task ${SLURM_ARRAY_TASK_ID}) - nothing to do, exiting cleanly."
    exit 0
fi

PAIR_FILE=$(echo "$MANIFEST_ROW" | cut -d',' -f1)
COMBINED_TRUE_LR_FILE=$(echo "$MANIFEST_ROW" | cut -d',' -f2)
UNRELATED_POOL_FILE=$(echo "$MANIFEST_ROW" | cut -d',' -f3)

echo "============================================================="
echo "  Pair Ranking Job (manifest-driven)"
echo "============================================================="
echo "  SLURM job ID          : ${SLURM_JOB_ID}"
echo "  Array task ID         : ${SLURM_ARRAY_TASK_ID}"
echo "  Manifest row (file)   : ${FILE_INDEX}"
echo "  Within-file chunk     : ${CHUNK_INDEX}"
echo "  Pair file             : ${PAIR_FILE}"
echo "  Combined true LR file : ${COMBINED_TRUE_LR_FILE}"
echo "  Unrelated pool file   : ${UNRELATED_POOL_FILE}"
echo "  Chunk size            : ${CHUNK_SIZE}"
echo "  Output dir            : ${OUTPUT_DIR}"
echo "  Started at            : $(date)"
echo "============================================================="

Rscript code/run_focal_ranking.R \
  "${PAIR_FILE}" \
  "${COMBINED_TRUE_LR_FILE}" \
  "${UNRELATED_POOL_FILE}" \
  "${CHUNK_INDEX}" \
  "${CHUNK_SIZE}" \
  "${OUTPUT_DIR}"

if [ $? -eq 0 ]; then
    echo "SUCCESS at $(date)"
else
    echo "ERROR at $(date)"
    exit 1
fi