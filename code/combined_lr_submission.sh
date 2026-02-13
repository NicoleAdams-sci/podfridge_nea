#!/bin/bash
# Combined LR Calculation Job Submission Script
# Simple cases - auto-generate file list from chunk range
	# sbatch code/combined_lr_submission.sh                    # All chunks, default list
	# sbatch code/combined_lr_submission.sh 1..20               # Chunks 1-20, default list
	# sbatch code/combined_lr_submission.sh 21..100             # Chunks 21-100, default list
# Advanced - provide custom file list (ignores chunk range for list generation)
	# sbatch code/combined_lr_submission.sh - my_custom_list.txt
	# sbatch code/combined_lr_submission.sh none my_custom_list.txt


# Use provided file list or create default in output directory
if [ "$1" = "-" ] || [ "$1" = "none" ]; then
    # Custom file list provided, use it as-is
    FILE_LIST=${2:-output/combined_lr_file_list.txt}
else
    # Auto-generate from chunk range
    CHUNK_RANGE=${1:-*}
    
    # Validate chunk range format - check for single dash (common mistake)
    if [[ "$CHUNK_RANGE" =~ ^[0-9]+-[0-9]+$ ]]; then
        echo "ERROR: Invalid chunk range format: $CHUNK_RANGE"
        echo "Use '..' instead of '-' for ranges (e.g., 1..20 not 1-20)"
        exit 1
    fi
    
    FILE_LIST="output/combined_lr_file_list.txt"
   
    if [ "$CHUNK_RANGE" = "*" ]; then
        # Default - match all chunks
        ls output/LR/LR_*_chunk*_*.csv > "$FILE_LIST" 2>/dev/null
    else
        # Specific range - use eval for brace expansion
        eval ls output/LR/LR_*_chunk{${CHUNK_RANGE}}_*.csv > "$FILE_LIST" 2>/dev/null
    fi
fi

# Count files and update SLURM script with file count AND file list path
FILE_COUNT=$(wc -l < "$FILE_LIST")
sed -e "s/FILE_COUNT_PLACEHOLDER/$FILE_COUNT/" \
    -e "s|FILE_LIST_PLACEHOLDER|$FILE_LIST|" \
    code/combined_lr.sh > code/combined_lr_ready.sh

# Submit job
sbatch code/combined_lr_ready.sh