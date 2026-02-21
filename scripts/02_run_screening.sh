#!/bin/bash

# --- RunPod Phase 1 Execution Script ---
# This script executes the high-throughput 'True Cas' test utilizing 
# the lightweight Protenix-Mini model to drastically save on compute costs.

JSON_DIR="../jsons"
# Output under phase1_screening so evolution_orchestrator can find PDBs
OUTPUT_DIR="../outputs/phase1_screening"
# Path to your localized databases to prevent IP leakage
DB_PATH="../databases" 

mkdir -p "$OUTPUT_DIR"

# Loop through all generated JSON configurations
for json_file in "$JSON_DIR"/*.json; do
    # Extract filename without extension for logging
    base_name=$(basename "$json_file" .json)
    echo "=================================================="
    echo "Processing Hit: $base_name"
    echo "=================================================="

    # Step 1: Strictly Localized Preprocessing
    # Generate the deep RNA MSAs needed for the crRNA and target RNA.
    # We use local paths to ensure proprietary sequences never hit an external API.
    echo "[1/2] Running localized MSA & Template Prep..."
    protenix inputprep \
        --input "$json_file" \
        --out_dir "$OUTPUT_DIR/${base_name}_prep" \
        --seqres_database_path "$DB_PATH" \
        > "$OUTPUT_DIR/${base_name}_prep.log" 2>&1

    # Check if prep succeeded
    if [ $? -ne 0 ]; then
        echo "Error: inputprep failed for $base_name. Check logs."
        continue
    fi

    # Step 2: Rapid Complex Assembly (Mini Model Screening)
    # Applying strict kernel optimizations and bfloat16 for the GPU cluster
    echo "[2/2] Running Protenix-Mini prediction..."
    protenix predict \
        --input "$OUTPUT_DIR/${base_name}_prep/$(basename "$json_file")" \
        --out_dir "$OUTPUT_DIR/${base_name}_pred" \
        --model_name "protenix_mini_default_v0.5.0" \
        --dtype bf16 \
        --enable_cache true \
        --enable_fusion true \
        --trimul_kernel true \
        > "$OUTPUT_DIR/${base_name}_pred.log" 2>&1

    echo "Completed $base_name. Outputs saved to $OUTPUT_DIR/${base_name}_pred/"
done

echo "Phase 1 High-Throughput Screening Complete."