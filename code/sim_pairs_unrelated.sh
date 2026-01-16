#!/bin/bash
#SBATCH --job-name=sim_unrelated
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
##SBATCH --mem=1G
#SBATCH --time=02:00:00
#SBATCH --array=1-400  # 5 populations × 80 chunks (was 450)
#SBATCH --output=logs/unrel_pairs_%A_%a.out
#SBATCH --error=logs/unrel_pairs_%A_%a.err

# Usage: sbatch code/sim_pairs_unrelated.sh

# Load required modules
module load Rtidyverse


POPULATIONS=("AfAm" "Cauc" "Hispanic" "Asian" "all")
CHUNKS_PER_POP=80  # (was 90)
N_PAIRS=1000

POP_INDEX=$(( (SLURM_ARRAY_TASK_ID - 1) / CHUNKS_PER_POP ))
CHUNK_INDEX=$(( (SLURM_ARRAY_TASK_ID - 1) % CHUNKS_PER_POP ))

POPULATION=${POPULATIONS[$POP_INDEX]}
#CHUNK_NUM=$(( CHUNK_INDEX + 11 ))  # Start at chunk 11 since you have 1-10
CHUNK_NUM=$(( CHUNK_INDEX + 21 ))  # Start at chunk 21 -- after two rounds of related pair runs

Rscript code/sim_pairs.R ${POPULATION} "unrelated" ${N_PAIRS} ${CHUNK_NUM}