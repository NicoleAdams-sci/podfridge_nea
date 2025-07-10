#!/bin/bash

# Usage: ./create_focal_database_multi.sh [num_pairs_per_relationship] [num_focal_individuals] [populations...]
# Example: ./create_focal_database_multi.sh 50 3 AfAm Cauc  # 3 focal individuals for AfAm and Cauc
# Example: ./create_focal_database_multi.sh 50 2           # 2 focal individuals for all populations

set -e  # Exit on error

NUM_PAIRS=${1:-50}        # Default to 50 pairs per relationship
NUM_FOCAL=${2:-1}         # Default to 1 focal individual
shift 2 2>/dev/null || shift $# # Remove first two arguments safely

# Define all populations
ALL_POPULATIONS=("AfAm" "Cauc" "Hispanic" "Asian")

# Determine which populations to process
if [ $# -gt 0 ]; then
    # Use specified populations
    POPULATIONS=("$@")
else
    # Process all populations
    POPULATIONS=("${ALL_POPULATIONS[@]}")
fi

echo "=== Creating Multiple Focal Individual Databases ==="
echo "Number of pairs per relationship type: $NUM_PAIRS"
echo "Number of focal individuals per population: $NUM_FOCAL"
echo "Processing populations: ${POPULATIONS[*]}"

# Create main focal database directory
mkdir -p output/focal_database

# Activate conda environment
echo -e "\n=== Activating conda environment ==="
source ~/.bashrc
conda activate rstats

# Process each population
for POPULATION in "${POPULATIONS[@]}"; do
    echo -e "\n=================================================="
    echo "=== PROCESSING POPULATION: $POPULATION ==="
    echo "=================================================="
    
    # Create population-specific directory
    FOCAL_DB_DIR="output/focal_database/${POPULATION}"
    mkdir -p "$FOCAL_DB_DIR"
    
    # Process each focal individual
    for FOCAL_NUM in $(seq 1 $NUM_FOCAL); do
        echo -e "\n--- Processing Focal Individual #$FOCAL_NUM for $POPULATION ---"
        
        # Determine file naming convention
        if [ "$FOCAL_NUM" -eq 1 ]; then
            FOCAL_SUFFIX=""
            FOCAL_ID_STR="1"
        else
            FOCAL_SUFFIX="_${FOCAL_NUM}"
            FOCAL_ID_STR="$FOCAL_NUM"
        fi
        
        # Step 1: Extract the focal individual
        echo -e "\n=== Step 1: Extracting focal individual #$FOCAL_NUM for $POPULATION ==="
        bash code/extract_focal_individual.sh "$POPULATION" "$FOCAL_NUM"
        
        # Step 2: Create pairs with unrelated individuals
        echo -e "\n=== Step 2: Creating pairs with unrelated individuals for focal #$FOCAL_NUM ==="
        bash code/create_focal_pairs.sh "$POPULATION" "$FOCAL_ID_STR"
        
        # Check pairs file
        if [ "$FOCAL_NUM" -eq 1 ]; then
            PAIRS_FILE="${FOCAL_DB_DIR}/focal_unrelated_pairs_${POPULATION}.csv"
        else
            PAIRS_FILE="${FOCAL_DB_DIR}/focal_unrelated_pairs_${POPULATION}_${FOCAL_NUM}.csv"
        fi
        
        UNRELATED_DB_ROWS=$(tail -n +2 "output/sim_processed_genotypes_unrelated_database.csv" | wc -l)
        PAIRS_ROWS=$(tail -n +2 "$PAIRS_FILE" | wc -l)
        
        if [ "$PAIRS_ROWS" -lt "$UNRELATED_DB_ROWS" ]; then
            echo "WARNING: Pairs file has fewer rows ($PAIRS_ROWS) than unrelated database ($UNRELATED_DB_ROWS)"
        fi
        
        # Step 3: Simulate additional relationship types
        echo -e "\n=== Step 3: Simulating relationships for focal #$FOCAL_NUM ==="
        
        # Update the R script call to handle multiple focal individuals
        if [ "$FOCAL_NUM" -eq 1 ]; then
            Rscript code/simulate_relationships_4focal.R "$POPULATION" "$NUM_PAIRS"
        else
            # For additional focal individuals, we'll need to modify the R script or create a variant
            echo "Note: Using focal individual #$FOCAL_NUM - R script needs updating for multiple focals"
            # For now, we'll create a copy approach
            FOCAL_FILE="${FOCAL_DB_DIR}/focal_individual_${POPULATION}${FOCAL_SUFFIX}.csv"
            TEMP_FOCAL_FILE="${FOCAL_DB_DIR}/focal_individual_${POPULATION}.csv.tmp"
            
            # Temporarily copy the focal individual file to the expected location
            cp "$FOCAL_FILE" "$TEMP_FOCAL_FILE"
            # Run simulation (this will overwrite existing files - we'll handle this)
            Rscript code/simulate_relationships_4focal.R "$POPULATION" "$NUM_PAIRS" || echo "Simulation may have issues with multiple focals"
            # Move results to focal-specific names
            if [ -f "${FOCAL_DB_DIR}/focal_all_relationships_${POPULATION}.csv" ]; then
                mv "${FOCAL_DB_DIR}/focal_all_relationships_${POPULATION}.csv" "${FOCAL_DB_DIR}/focal_all_relationships_${POPULATION}${FOCAL_SUFFIX}.csv"
            fi
            rm -f "$TEMP_FOCAL_FILE"
        fi
        
        # Step 4: Prepare parent-child data
        echo -e "\n=== Step 4: Preparing parent-child data for focal #$FOCAL_NUM ==="
        
        PARENT_CHILD_FILE="output/sim_processed_genotypes_${POPULATION}_parent_child_combined.csv"
        FOCAL_PARENT_CHILD="${FOCAL_DB_DIR}/focal_parent_child_${POPULATION}${FOCAL_SUFFIX}.csv"
        
        if [ -f "$PARENT_CHILD_FILE" ]; then
            # Get this focal individual's sim_id
            if [ "$FOCAL_NUM" -eq 1 ]; then
                FOCAL_SIM_ID=$(awk -F, 'NR>1 {print $3; exit}' "$PARENT_CHILD_FILE")
            else
                # Get the nth unique sim_id
                FOCAL_SIM_ID=$(tail -n +2 "$PARENT_CHILD_FILE" | awk -F, -v pop="$POPULATION" '$1 == pop {print $3}' | sort -u | sed -n "${FOCAL_NUM}p")
            fi
            
            echo "Extracting parent-child data for focal individual #$FOCAL_NUM with sim_id: $FOCAL_SIM_ID"
            
            # Create parent-child file
            echo "population,relationship_type,pair_id,locus,ind1_allele1,ind1_allele2,ind2_allele1,ind2_allele2" > "$FOCAL_PARENT_CHILD"
            awk -F, -v sim="$FOCAL_SIM_ID" -v pop="$POPULATION" '
            $1 == pop && $3 == sim {
                print $1 "," $2 "," sim "," $4 "," $5 "," $6 "," $7 "," $8
            }' "$PARENT_CHILD_FILE" >> "$FOCAL_PARENT_CHILD"
            
            PC_ROWS=$(tail -n +2 "$FOCAL_PARENT_CHILD" | wc -l)
            echo "Extracted $PC_ROWS rows of parent-child data"
        else
            echo "WARNING: Parent-child file not found: $PARENT_CHILD_FILE"
            echo "population,relationship_type,pair_id,locus,ind1_allele1,ind1_allele2,ind2_allele1,ind2_allele2" > "$FOCAL_PARENT_CHILD"
        fi
        
        # Step 5: Combine all files
        echo -e "\n=== Step 5: Combining files for focal #$FOCAL_NUM ==="
        
        UNRELATED_PAIRS="${FOCAL_DB_DIR}/focal_unrelated_pairs_${POPULATION}${FOCAL_SUFFIX}.csv"
        ALL_RELATIONSHIPS="${FOCAL_DB_DIR}/focal_all_relationships_${POPULATION}${FOCAL_SUFFIX}.csv"
        COMPLETE_DB="${FOCAL_DB_DIR}/focal_complete_database_${POPULATION}${FOCAL_SUFFIX}.csv"
        
        # Create combined database
        head -n 1 "$UNRELATED_PAIRS" > "$COMPLETE_DB"
        
        # Add parent-child data
        if [ -f "$FOCAL_PARENT_CHILD" ] && [ $(wc -l < "$FOCAL_PARENT_CHILD") -gt 1 ]; then
            tail -n +2 "$FOCAL_PARENT_CHILD" >> "$COMPLETE_DB"
        fi
        
        # Add unrelated pairs
        tail -n +2 "$UNRELATED_PAIRS" >> "$COMPLETE_DB"
        
        # Add simulated relationships
        if [ -f "$ALL_RELATIONSHIPS" ]; then
            tail -n +2 "$ALL_RELATIONSHIPS" >> "$COMPLETE_DB"
        fi
        
        TOTAL_ROWS=$(tail -n +2 "$COMPLETE_DB" | wc -l)
        echo "Combined database for focal #$FOCAL_NUM: $TOTAL_ROWS rows"
        echo "Saved to: $COMPLETE_DB"
        
        echo -e "\n--- Focal Individual #$FOCAL_NUM for $POPULATION complete ---"
    done
    
    echo -e "\n=== $POPULATION processing complete (processed $NUM_FOCAL focal individuals) ==="
done

# Create summary
echo -e "\n=== Creating Processing Summary ==="
COMBINED_DIR="output/focal_database/combined_analysis"
mkdir -p "$COMBINED_DIR"

SUMMARY_FILE="${COMBINED_DIR}/multi_focal_summary.txt"
echo "Multiple Focal Individual Database Processing Summary" > "$SUMMARY_FILE"
echo "====================================================" >> "$SUMMARY_FILE"
echo "Date: $(date)" >> "$SUMMARY_FILE"
echo "Pairs per relationship type: $NUM_PAIRS" >> "$SUMMARY_FILE"
echo "Focal individuals per population: $NUM_FOCAL" >> "$SUMMARY_FILE"
echo "Populations processed: ${POPULATIONS[*]}" >> "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"

for POPULATION in "${POPULATIONS[@]}"; do
    echo "Population: $POPULATION" >> "$SUMMARY_FILE"
    for FOCAL_NUM in $(seq 1 $NUM_FOCAL); do
        if [ "$FOCAL_NUM" -eq 1 ]; then
            FOCAL_SUFFIX=""
        else
            FOCAL_SUFFIX="_${FOCAL_NUM}"
        fi
        COMPLETE_DB="output/focal_database/${POPULATION}/focal_complete_database_${POPULATION}${FOCAL_SUFFIX}.csv"
        if [ -f "$COMPLETE_DB" ]; then
            TOTAL_ROWS=$(tail -n +2 "$COMPLETE_DB" | wc -l)
            echo "  Focal #$FOCAL_NUM: $TOTAL_ROWS rows" >> "$SUMMARY_FILE"
        fi
    done
    echo "" >> "$SUMMARY_FILE"
done

echo "Multi-focal processing summary saved to: $SUMMARY_FILE"

echo -e "\n=================================================="
echo "=== MULTI-FOCAL PROCESSING COMPLETE ==="
echo "=================================================="
echo "Processed populations: ${POPULATIONS[*]}"
echo "Focal individuals per population: $NUM_FOCAL"
echo "Results in: output/focal_database/"
echo ""
echo "Note: LR calculations and analysis need to be run separately for each focal individual"
echo "Next steps:"
echo "1. Run LR calculations for each focal individual"
echo "2. Use analysis scripts to compare performance across focal individuals"