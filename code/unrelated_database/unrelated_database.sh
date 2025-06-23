#!/bin/bash

### make "database" of unrelated individuals pulled from the simulated genotypes
### across populations at specified population proportions


# Extract the total sim_ids parameter (each sim_id will have 29 loci)
TOTAL_SIM_IDS=$1
shift  # Remove the first parameter (total rows), leaving only directories


# Populations to process
POPULATIONS=("AfAm" "Cauc" "Hispanic" "Asian")

# Predefined sampling numbers for unrelated rows
UNRELATED_SAMPLING=(
    "AfAm:20"
    "Cauc:50"
    "Asian:10"
    "Hispanic:20"
)

# Ensure we're in the correct output directory
if [ ! -d "output/" ]; then
    echo "Error: Cannot find ../output directory"
    exit 1
fi
cd output/ || exit 1


# Predefined header
FULL_HEADER="population,relationship_type,sim_id,locus,ind1_allele1,ind1_allele2,ind2_allele1,ind2_allele2,shared_alleles,genotype_match,LR_AfAm,LR_Cauc,LR_Hispanic,LR_Asian"
DATABASE_HEADER="population,relationship_type,db_id,locus,ind2_allele1,ind2_allele2"

# Sample unrelated rows for "database"

# Create/clear the combined unrelated sample file with new name
> sim_processed_genotypes_unrelated_database.csv

# Write header to the sampled unrelated file
echo "$DATABASE_HEADER" > sim_processed_genotypes_unrelated_database.csv


# Process sampling for each population
for SAMPLE in "${UNRELATED_SAMPLING[@]}"; do
    # Split population and percentage
    POP=${SAMPLE%%:*}
    PERCENT=${SAMPLE##*:}
    
    # Calculate number of sim_ids to sample based on percentage
    NUM_SIM_IDS=$(( (TOTAL_SIM_IDS * PERCENT) / 100 ))
    
    # Input unrelated file
    UNRELATED_FILE="sim_processed_genotypes_${POP}_unrelated_combined.csv"
    
    # Create a temporary file to store unique sim_ids
    TMP_SIM_IDS=$(mktemp)
    
    # Extract unique sim_ids from the unrelated file (skip header)
    tail -n +2 "$UNRELATED_FILE" | awk -F, '{print $3}' | sort -u > "$TMP_SIM_IDS"
    
    # Count total unique sim_ids
    TOTAL_UNIQUE_SIMS=$(wc -l < "$TMP_SIM_IDS")
    echo "Population $POP: Sampling $NUM_SIM_IDS sim_ids from $TOTAL_UNIQUE_SIMS available"
    
    # Sample the required number of sim_ids
    SAMPLED_SIM_IDS=$(mktemp)
    shuf -n "$NUM_SIM_IDS" "$TMP_SIM_IDS" > "$SAMPLED_SIM_IDS"
    
    # For each sampled sim_id, extract all corresponding rows (using exact matching)
    while read -r SIM_ID; do
        # Use awk for exact field matching instead of grep
        #awk -F, -v sid="$SIM_ID" -v pop="$POP" '$1 == pop && $3 == sid' "$UNRELATED_FILE" >> sim_processed_genotypes_unrelated_database.csv # keep full row
        awk -F, -v sid="$SIM_ID" -v pop="$POP" '$1 == pop && $3 == sid {print $1","$2","$3","$4","$5","$6}' "$UNRELATED_FILE" >> sim_processed_genotypes_unrelated_database.csv # only keep indiv1
    done < "$SAMPLED_SIM_IDS"
    
    # Clean up temporary files
    rm "$TMP_SIM_IDS" "$SAMPLED_SIM_IDS"
done

# Count actual number of rows and sim_ids in the final file (excluding header)
ACTUAL_ROWS=$(( $(wc -l < sim_processed_genotypes_unrelated_database.csv) - 1 ))
ACTUAL_SIM_IDS=$(tail -n +2 sim_processed_genotypes_unrelated_database.csv | awk -F, '{print $1 "," $3}' | sort -u | wc -l)
EXPECTED_ROWS=$(( TOTAL_SIM_IDS * 29 ))
# Check if we got the expected number of rows
if [ "$ACTUAL_ROWS" -ne "$EXPECTED_ROWS" ]; then
    echo "WARNING: Expected $EXPECTED_ROWS rows (${TOTAL_SIM_IDS} sim_ids × 29 loci), but got $ACTUAL_ROWS rows"
    echo "Analyzing rows per sim_id:"
    
    # Create a temporary analysis file
    ANALYSIS_FILE=$(mktemp)
    tail -n +2 sim_processed_genotypes_unrelated_database.csv | awk -F, '{print $1 "," $3}' | sort | uniq -c | sort -nr > "$ANALYSIS_FILE"
    
    # Show top 10 sim_ids with counts
    echo "Top sim_ids by row count:"
    head -10 "$ANALYSIS_FILE"
    
    # Check if any sim_ids have more or less than 29 loci
    echo "Sim_ids with incorrect number of loci:"
    awk '$1 != 29 {print $0}' "$ANALYSIS_FILE" | head -10
    
    # Clean up
    rm "$ANALYSIS_FILE"
fi

echo "Unrelated database created successfully:"
echo "- Requested sim_ids: $TOTAL_SIM_IDS"
echo "- Actual sim_ids: $ACTUAL_SIM_IDS"
echo "- Expected rows: $EXPECTED_ROWS"
echo "- Actual rows: $ACTUAL_ROWS"
echo "- Output file: sim_processed_genotypes_unrelated_database.csv"