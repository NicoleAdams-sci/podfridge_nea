#!/bin/bash

# Populations to process
POPULATIONS=("AfAm" "Cauc" "Hispanic" "Asian")

# Relationship types to extract
RELATIONSHIP_TYPES=("parent_child" "full_siblings" "half_siblings" "cousins" "second_cousins" "unrelated")

# Predefined header
FULL_HEADER="population,relationship_type,sim_id,locus,ind1_allele1,ind1_allele2,ind2_allele1,ind2_allele2,shared_alleles,genotype_match,LR_AfAm,LR_Cauc,LR_Hispanic,LR_Asian"

# Ensure we're in the correct output directory
if [ ! -d "output/" ]; then
    echo "Error: Cannot find ../output directory"
    exit 1
fi
cd output/ || exit 1

# Expand directories explicitly to handle wildcards
DIRS=($@)

# Process populations for each directory
for DIR in "${DIRS[@]}"; do
    # Ensure directory exists and is accessible
    if [ ! -d "$DIR" ]; then
        echo "Error: Cannot access directory $DIR"
        continue
    fi
    
    # Change to directory and process
    cd "$DIR" || continue
    
    # Extract population and relationship-specific data
    for POP in "${POPULATIONS[@]}"; do
        for REL in "${RELATIONSHIP_TYPES[@]}"; do
            grep "$POP" sim_processed_genotypes_*.csv | grep "$REL" > "${POP}_${REL}_processed_genotypes.csv" 2>/dev/null
        done
    done
    
    # Return to parent directory
    cd ..
done

# Combine population and relationship files across directories
for POP in "${POPULATIONS[@]}"; do
    for REL in "${RELATIONSHIP_TYPES[@]}"; do
        # Create/clear combined files
        > "sim_processed_genotypes_${POP}_${REL}_combined.csv"
        
        # Write predefined header
        echo "$FULL_HEADER" > "sim_processed_genotypes_${POP}_${REL}_combined.csv"
        
        # Append data from all directories
        for DIR in "${DIRS[@]}"; do
            POPULATION_REL_FILE="$DIR/${POP}_${REL}_processed_genotypes.csv"
            if [ -f "$POPULATION_REL_FILE" ]; then
                # Append all lines (no header skipping needed)
                cat "$POPULATION_REL_FILE" >> "sim_processed_genotypes_${POP}_${REL}_combined.csv"
                
                # Remove individual population-relationship file after combining
                rm "$POPULATION_REL_FILE"
            fi
        done
    done
done

echo "Population and relationship data combined and temporary files removed successfully"
