#!/bin/bash

# Usage: ./create_focal_pairs.sh <population>
# Example: ./create_focal_pairs.sh AfAm

POPULATION=$1

if [ -z "$POPULATION" ]; then
    echo "Error: Population parameter is required"
    echo "Usage: ./create_focal_pairs.sh <population>"
    exit 1
fi

# Define paths
UNRELATED_DB="output/sim_processed_genotypes_unrelated_database.csv"
FOCAL_FILE="output/focal_database/${POPULATION}/focal_individual_${POPULATION}.csv"
OUTPUT_FILE="output/focal_database/${POPULATION}/focal_unrelated_pairs_${POPULATION}.csv"

# Check if the required files exist
if [ ! -f "$UNRELATED_DB" ]; then
    echo "Error: Cannot find unrelated database: $UNRELATED_DB"
    exit 1
fi

if [ ! -f "$FOCAL_FILE" ]; then
    echo "Error: Cannot find focal individual file: $FOCAL_FILE"
    exit 1
fi

# Create output directory if needed
mkdir -p $(dirname "$OUTPUT_FILE")

# Define headers
PAIR_HEADER="population,relationship_type,pair_id,locus,focal_allele1,focal_allele2,ind2_allele1,ind2_allele2"

# Write header to output file
echo "$PAIR_HEADER" > "$OUTPUT_FILE"

# Count total rows in the unrelated database
TOTAL_UNRELATED_ROWS=$(tail -n +2 "$UNRELATED_DB" | wc -l)
echo "Total rows in unrelated database: $TOTAL_UNRELATED_ROWS"

# Get the focal individual's data into an associative array
declare -A focal_alleles
while IFS=, read -r pop focal_id locus allele1 allele2; do
    # Skip header
    if [ "$pop" = "population" ]; then continue; fi
    
    # Store alleles by locus
    focal_alleles["$locus"]="$allele1,$allele2"
done < "$FOCAL_FILE"

echo "Loaded focal individual data for ${#focal_alleles[@]} loci"

# Process all unrelated individuals regardless of population
pair_id=1
while IFS=, read -r pop rel_type sim_id locus ind1_a1 ind1_a2; do
    # Skip header
    if [ "$pop" = "population" ]; then continue; fi
    
    # Check if we have focal alleles for this locus
    if [ -n "${focal_alleles[$locus]}" ]; then
        IFS=, read -r focal_a1 focal_a2 <<< "${focal_alleles[$locus]}"
        
        # Get pair_id from sim_id to maintain 1:1 mapping
        # Extract the original sim_id number to use as pair_id
        
        # Write pair data: focal first, then unrelated individual
        echo "$pop,unrelated_focal,$sim_id,$locus,$focal_a1,$focal_a2,$ind1_a1,$ind1_a2" >> "$OUTPUT_FILE"
    else
        echo "Warning: No focal allele data for locus $locus"
    fi
done < "$UNRELATED_DB"

# Count the rows in the output file
ROWS=$(tail -n +2 "$OUTPUT_FILE" | wc -l)

echo "Created pairs between unrelated individuals and the focal individual"
echo "Output contains $ROWS rows (expected to match the unrelated database)"
echo "Output saved to: $OUTPUT_FILE"

# Check if we have fewer rows than expected
if [ "$ROWS" -lt "$TOTAL_UNRELATED_ROWS" ]; then
    echo "WARNING: Output has fewer rows ($ROWS) than the unrelated database ($TOTAL_UNRELATED_ROWS)"
    echo "This could be because some loci in the unrelated database are not present in the focal individual data"
    
    # Find which loci are missing
    echo "Checking which loci might be missing..."
    MISSING_LOCI=$(mktemp)
    
    # Extract unique loci from unrelated database
    tail -n +2 "$UNRELATED_DB" | awk -F, '{print $4}' | sort -u > "$MISSING_LOCI.unrelated"
    
    # Extract loci that have focal alleles
    for key in "${!focal_alleles[@]}"; do
        echo "$key" >> "$MISSING_LOCI.focal"
    done
    sort "$MISSING_LOCI.focal" > "$MISSING_LOCI.focal.sorted"
    
    # Find loci in unrelated but not in focal
    comm -23 "$MISSING_LOCI.unrelated" "$MISSING_LOCI.focal.sorted" > "$MISSING_LOCI.diff"
    
    echo "Loci in unrelated database but missing in focal individual:"
    cat "$MISSING_LOCI.diff"
    
    # Clean up
    rm "$MISSING_LOCI" "$MISSING_LOCI.unrelated" "$MISSING_LOCI.focal" "$MISSING_LOCI.focal.sorted" "$MISSING_LOCI.diff"
fi
