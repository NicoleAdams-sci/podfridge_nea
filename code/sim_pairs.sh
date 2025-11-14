#!/bin/bash
#SBATCH --job-name=sim_pairs
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
##SBATCH --mem=1G
#SBATCH --time=00:30:00
#SBATCH --array=1-300   # 5 populations × 6 relationships × 10 chunks
#SBATCH --output=logs/sim_pairs_%A_%a.out
#SBATCH --error=logs/sim_pairs_%A_%a.err

# Usage: sbatch code/sim_pairs.sh 100 [100 pairs of each pop x relationship]
# Usage: sbatch  --wrap="module load Rtidyverse; Rscript code/sim_pairs_test.R Cauc full_siblings 100" [100 pairs of just Cauc full_siblings]

# Create logs directory if it doesn't exist
mkdir -p logs

# Get N_PAIRS parameter with default of 1
N_PAIRS=${1:-1}

# Load required modules
module load Rtidyverse

# Define populations and relationships
POPULATIONS=("AfAm" "Cauc" "Hispanic" "Asian" "all")
RELATIONSHIPS=("parent_child" "full_siblings" "half_siblings" "cousins" "second_cousins" "unrelated")

# Calculate which combination and chunk based on SLURM_ARRAY_TASK_ID
# For 10 chunks of 1000 pairs each per combination
CHUNKS_PER_COMBO=10
N_PAIRS=1000

# Calculate combination index (0-23) and chunk index (0-9)
COMBO_INDEX=$(( (SLURM_ARRAY_TASK_ID - 1) / CHUNKS_PER_COMBO ))
CHUNK_INDEX=$(( (SLURM_ARRAY_TASK_ID - 1) % CHUNKS_PER_COMBO ))

# Calculate population and relationship from combination index
POP_INDEX=$(( COMBO_INDEX / 6 ))
REL_INDEX=$(( COMBO_INDEX % 6 ))

POPULATION=${POPULATIONS[$POP_INDEX]}
RELATIONSHIP=${RELATIONSHIPS[$REL_INDEX]}
CHUNK_NUM=$(( CHUNK_INDEX + 1 ))  # Make it 1-based for display


echo "Running job ${SLURM_ARRAY_TASK_ID}: Population=${POPULATION}, Relationship=${RELATIONSHIP}, Chunk=${CHUNK_NUM}"
echo "Started at: $(date)"

# Run the R simulate pairs (module 3) script
Rscript code/sim_pairs.R ${POPULATION} ${RELATIONSHIP} ${N_PAIRS} ${CHUNK_NUM}


echo "Completed at: $(date)"