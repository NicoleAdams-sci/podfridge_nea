#!/bin/bash

# Usage: ./extract_focal_individual.sh <population> [focal_individual_number]
# Example: ./extract_focal_individual.sh AfAm
# Example: ./extract_focal_individual.sh AfAm 5

POPULATION=$1
FOCAL_NUMBER=${2:-1}  # Default to first individual

if [ -z "$POPULATION" ]; then
    echo "Error: Population parameter is required"
    echo "Usage: ./extract_focal_individual.sh <population> [focal_individual_number]"
    echo "  focal_individual_number: Which individual to extract (1=first, 2=second, etc.)"
    exit 1
fi

# Define paths with new directory structure
PARENT_CHILD_FILE="output/sim_processed_genotypes_${POPULATION}_parent_child_combined.csv"
FOCAL_DB_DIR="output/focal_database/${POPULATION}"

# Support multiple focal individuals
if [ "$FOCAL_NUMBER" -eq 1 ]; then
    OUTPUT_FILE="${FOCAL_DB_DIR}/focal_individual_${POPULATION}.csv"
else
    OUTPUT_FILE="${FOCAL_DB_DIR}/focal_individual_${POPULATION}_${FOCAL_NUMBER}.csv"
fi

# Check if the parent-child file exists
if [ ! -f "$PARENT_CHILD_FILE" ]; then
    echo "Error: Cannot find parent-child file: $PARENT_CHILD_FILE"
    exit 1
fi

# Create the focal database directory structure
mkdir -p "$FOCAL_DB_DIR"

# Define headers
ORIGINAL_HEADER=$(head -n 1 "$PARENT_CHILD_FILE")
FOCAL_HEADER="population,focal_id,locus,focal_allele1,focal_allele2"

# Write header to output file
echo "$FOCAL_HEADER" > "$OUTPUT_FILE"

# Get all unique sim_ids from the parent-child file
TEMP_SIM_IDS=$(mktemp)
tail -n +2 "$PARENT_CHILD_FILE" | awk -F, -v pop="$POPULATION" '$1 == pop {print $3}' | sort -u > "$TEMP_SIM_IDS"

TOTAL_AVAILABLE=$(wc -l < "$TEMP_SIM_IDS")

echo "Found $TOTAL_AVAILABLE unique individuals in parent-child data for $POPULATION"

if [ "$FOCAL_NUMBER" -gt "$TOTAL_AVAILABLE" ]; then
    echo "Error: Requested focal individual number $FOCAL_NUMBER, but only $TOTAL_AVAILABLE individuals available"
    rm "$TEMP_SIM_IDS"
    exit 1
fi

# Extract the specified focal individual (nth line from the sorted list)
FOCAL_SIM_ID=$(sed -n "${FOCAL_NUMBER}p" "$TEMP_SIM_IDS")
rm "$TEMP_SIM_IDS"

echo "Selected focal individual #$FOCAL_NUMBER with sim_id: $FOCAL_SIM_ID"

# Extract all loci for this focal individual (ind1), transforming the data format
awk -F, -v sim="$FOCAL_SIM_ID" -v pop="$POPULATION" '
$1 == pop && $3 == sim {
    print $1 "," sim "," $4 "," $5 "," $6
}' "$PARENT_CHILD_FILE" >> "$OUTPUT_FILE"

# Count rows in the output file
ROWS=$(wc -l < "$OUTPUT_FILE")
EXPECTED_ROWS=30  # Header + 29 loci

if [ "$ROWS" -ne "$EXPECTED_ROWS" ]; then
    echo "Warning: Expected $EXPECTED_ROWS rows (1 header + 29 loci), but got $ROWS rows"
else
    echo "Successfully extracted focal individual with $((ROWS-1)) loci"
fi

echo "Focal individual data saved to: $OUTPUT_FILE"
