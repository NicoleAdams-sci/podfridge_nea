#!/bin/bash
#SBATCH --job-name=focal_family
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/focal_family_%j.out
#SBATCH --error=logs/focal_family_%j.err

# ==============================================================================
# Generate Focal Families — Shell Wrapper
# ==============================================================================
#
# Generates focal individuals and their relatives for a single population,
# saved to output/focal_ranking_test/.
#
# Usage — predefined family structure:
#   bash code/generate_focal_family.sh 10 Asian minimal
#   bash code/generate_focal_family.sh 10 Asian extended
#
# Usage — custom relationship counts (FAMILY_TYPE = "custom"):
#   bash code/generate_focal_family.sh 10 Asian custom "second_cousins=1"
#   bash code/generate_focal_family.sh 10 Asian custom "parent_child=1,full_siblings=1,cousins=1"
#   bash code/generate_focal_family.sh 10 Asian custom "half_siblings=1,second_cousins=1"
#
# Arguments:
#   $1 : N_FOCAL        — number of focal individuals (default: 10)
#   $2 : POPULATION     — focal individual population (default: Asian)
#   $3 : FAMILY_TYPE    — predefined type OR "custom" (default: minimal)
#                         predefined options: minimal, nuclear, large_family,
#                                            extended, siblings_only, parent_only
#   $4 : CUSTOM_COUNTS  — comma-separated relationship=count pairs (only used
#                         when FAMILY_TYPE="custom")
#                         e.g. "second_cousins=1" or "parent_child=1,cousins=2"
#
# Output:
#   output/focal_ranking_test/<pop>_focal_<structure>_<datetime>.csv
# ==============================================================================

mkdir -p logs
mkdir -p output/focal_ranking_test

module load Rtidyverse

N_FOCAL=${1:-10}
POPULATION=${2:-Asian}
FAMILY_TYPE=${3:-minimal}
CUSTOM_COUNTS=${4:-""}

echo "============================================================="
echo "  Generate Focal Families"
echo "============================================================="
echo "  Population      : ${POPULATION}"
echo "  N focal         : ${N_FOCAL}"
echo "  Family type     : ${FAMILY_TYPE}"
if [ "${FAMILY_TYPE}" = "custom" ]; then
    echo "  Custom counts   : ${CUSTOM_COUNTS}"
fi
echo "  Output directory: output/focal_ranking_test/"
echo "  Started at      : $(date)"
echo "============================================================="

Rscript code/generate_focal_family.R \
    ${N_FOCAL} \
    ${POPULATION} \
    ${FAMILY_TYPE} \
    "${CUSTOM_COUNTS}"

if [ $? -eq 0 ]; then
    echo "SUCCESS: Focal families generated at $(date)"
else
    echo "ERROR: Focal family generation failed at $(date)"
    exit 1
fi
