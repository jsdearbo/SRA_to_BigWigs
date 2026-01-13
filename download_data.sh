#!/bin/bash
set -euo pipefail

# ============================================================
# Configuration
# ============================================================
METADATA_FILE="./metadata.csv"
OUTPUT_DIR="/mnt/netfiles/jsdearbo/public_data/pipeline_v2/sra_files"

# Ensure tools are available
command -v prefetch >/dev/null 2>&1 || { echo "Error: 'prefetch' not found."; exit 1; }

# ============================================================
# Main Loop
# ============================================================

mkdir -p "$OUTPUT_DIR"

# Skip the header line (NR>1) and read comma-separated values
awk -F',' 'NR>1 {print $1, $2, $4}' "$METADATA_FILE" | while read -r cell_type run_id description; do
    
    # Create organization directory (optional, but good for structure)
    # We download everything to one SRA folder to keep prefetch happy, 
    # but we log what we are doing.
    
    echo "---------------------------------------------------"
    echo "Target: $cell_type | ID: $run_id"
    echo "Desc:   $description"
    
    target_file="$OUTPUT_DIR/$run_id/$run_id.sra"
    
    if [[ -f "$target_file" ]]; then
        echo "[SKIP] $run_id already exists."
    else
        echo "[DOWNLOADING] Fetching $run_id..."
        prefetch --max-size 50G -O "$OUTPUT_DIR" "$run_id"
    fi

done

echo "---------------------------------------------------"
echo "All downloads processed based on $METADATA_FILE."