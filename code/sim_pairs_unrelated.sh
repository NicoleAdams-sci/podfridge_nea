#!/bin/bash
#SBATCH --job-name=sim_unrelated
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
##SBATCH --mem=1G
#SBATCH --time=02:00:00
#SBATCH --array=1-400  # 5 populations × 80 chunks
#SBATCH --output=logs/unrel_pairs_%A_%a.out
#SBATCH --error=logs/unrel_pairs_%A_%a.err

# Usage: sbatch code/sim_pairs_unrelated.sh           # Default: 1000 pairs/chunk (80k total)
# Usage: sbatch code/sim_pairs_unrelated.sh 100       # Testing: 100 pairs/chunk (8k total)
# Usage: sbatch code/sim_pairs_unrelated.sh 500       # 500 pairs/chunk (40k total)

# Load required modules
module load Rtidyverse


POPULATIONS=("AfAm" "Cauc" "Hispanic" "Asian" "all")
CHUNKS_PER_POP=80
N_PAIRS=${1:-1000}  # Default to 1000 if no argument provided

POP_INDEX=$(( (SLURM_ARRAY_TASK_ID - 1) / CHUNKS_PER_POP ))
CHUNK_INDEX=$(( (SLURM_ARRAY_TASK_ID - 1) % CHUNKS_PER_POP ))

POPULATION=${POPULATIONS[$POP_INDEX]}
CHUNK_NUM=$(( CHUNK_INDEX + 21 ))  # Start at chunk 21 (after related pairs in chunks 1-20)

Rscript code/sim_pairs.R ${POPULATION} "unrelated" ${N_PAIRS} ${CHUNK_NUM}