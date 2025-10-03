#!/bin/bash
#SBATCH --job-name=sim_pairs
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
##SBATCH --mem=1G
#SBATCH --time=00:10:00
#SBATCH --array=1-24
#SBATCH --output=output/logs/sim_pairs_%A_%a.out
#SBATCH --error=output/logs/sim_pairs_%A_%a.err

# Usage: sbatch code/sim_pairs.sh 100 [100 pairs of each pop x relationship]
# Usage: sbatch  --wrap="module load Rtidyverse; Rscript code/sim_pairs_test.R Cauc full_siblings 100" [100 pairs of just Cauc full_siblings]

# Create logs directory if it doesn't exist
mkdir -p output/logs

# Get N_PAIRS parameter with default of 1
N_PAIRS=${1:-1}

# Load required modules
module load Rtidyverse

# Define populations and relationships
POPULATIONS=("AfAm" "Cauc" "Hispanic" "Asian")
RELATIONSHIPS=("parent_child" "full_siblings" "half_siblings" "cousins" "second_cousins" "unrelated")

# Calculate which combination to run based on SLURM_ARRAY_TASK_ID
# Total combinations: 4 populations × 6 relationships = 24 jobs
POP_INDEX=$(( (SLURM_ARRAY_TASK_ID - 1) / 6 ))
REL_INDEX=$(( (SLURM_ARRAY_TASK_ID - 1) % 6 ))

POPULATION=${POPULATIONS[$POP_INDEX]}
RELATIONSHIP=${RELATIONSHIPS[$REL_INDEX]}

echo "Running job ${SLURM_ARRAY_TASK_ID}: Population=${POPULATION}, Relationship=${RELATIONSHIP}"
echo "Started at: $(date)"

# Run the R simulate pairs (module 3) script
Rscript code/sim_pairs_test.R ${POPULATION} ${RELATIONSHIP} $N_PAIRS


echo "Completed at: $(date)"