#!/bin/bash

# --- RunPod Phase 1 Execution Script ---
# This script executes the high-throughput 'True Cas' test utilizing
# the lightweight Protenix-Mini model to drastically save on compute costs.
# Uses protenix msa (when available) for better quality; falls back to raw JSON.
#
# Processes all baselines. Validation filtering happens before evolution only.

log_ts() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

JSON_DIR="../jsons"
OUTPUT_DIR="../outputs/phase1_screening"
SKIP_MSA="${SKIP_MSA:-1}"

mkdir -p "$OUTPUT_DIR"

# Loop through all generated JSON configurations
for json_file in "$JSON_DIR"/*.json; do
    [ -f "$json_file" ] || continue
    base_name=$(basename "$json_file" .json)
    echo "=================================================="
    echo "Processing Hit: $base_name"
    echo "=================================================="

    # Step 1: Optional MSA (skipped by default - remote server is slow)
    PREDICT_INPUT="$json_file"
    USE_MSA="false"
    if [ "$SKIP_MSA" != "1" ]; then
        MSA_DIR="$OUTPUT_DIR/${base_name}_msa"
        mkdir -p "$MSA_DIR"
        echo "[1/2] Running MSA search (protenix msa)..."
        if protenix msa --input "$json_file" --out_dir "$MSA_DIR" > "$OUTPUT_DIR/${base_name}_msa.log" 2>&1; then
            MSA_JSON="$MSA_DIR/$(basename "$json_file")"
            if [ -f "$MSA_JSON" ]; then
                PREDICT_INPUT="$MSA_JSON"
                USE_MSA="true"
            fi
        else
            log_ts "MSA failed for $base_name, using raw JSON..."
        fi
    else
        log_ts "[1/2] Skipping MSA (SKIP_MSA=1), using raw JSON..."
    fi

    # Step 2: Protenix-Mini prediction
    log_ts "[2/2] Running Protenix-Mini prediction (may take 2-10 min per hit)..."
    protenix predict \
        --input "$PREDICT_INPUT" \
        --out_dir "$OUTPUT_DIR/${base_name}_pred" \
        --model_name "protenix_mini_default_v0.5.0" \
        --use_msa "$USE_MSA" \
        --use_default_params true \
        > "$OUTPUT_DIR/${base_name}_pred.log" 2>&1

    log_ts "Completed $base_name. Outputs saved to $OUTPUT_DIR/${base_name}_pred/"
done

log_ts "Phase 1 High-Throughput Screening Complete."