#!/bin/bash

# --- RunPod Phase 1 Execution Script ---
# This script executes the high-throughput 'True Cas' test utilizing
# the lightweight Protenix-Mini model to drastically save on compute costs.
# Uses protenix msa (when available) for better quality; falls back to raw JSON.
#
# Multi-GPU: set NUM_GPUS to parallelize across GPUs (e.g. NUM_GPUS=4).
# Skips baselines that already have outputs in OUTPUT_DIR.
# Processes all baselines. Validation filtering happens before evolution only.

log_ts() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

JSON_DIR="../jsons"
OUTPUT_DIR="../outputs/phase1_screening"
SKIP_MSA="${SKIP_MSA:-1}"
NUM_GPUS="${NUM_GPUS:-1}"

mkdir -p "$OUTPUT_DIR"

run_single() {
    local json_file="$1"
    local gpu_id="$2"
    local base_name
    base_name=$(basename "$json_file" .json)

    # Skip if outputs already exist (resume support)
    local pred_dir="$OUTPUT_DIR/${base_name}_pred"
    if [ -d "$pred_dir" ] && find "$pred_dir" -name '*.cif' -o -name '*.pdb' 2>/dev/null | head -1 | grep -q .; then
        log_ts "Skipping $base_name (outputs already exist)"
        return 0
    fi

    echo "=================================================="
    echo "[GPU $gpu_id] Processing Hit: $base_name"
    echo "=================================================="

    # Step 1: Optional MSA (skipped by default - remote server is slow)
    PREDICT_INPUT="$json_file"
    USE_MSA="false"
    if [ "$SKIP_MSA" != "1" ]; then
        MSA_DIR="$OUTPUT_DIR/${base_name}_msa"
        mkdir -p "$MSA_DIR"
        echo "[1/2] Running MSA search (protenix msa)..."
        if CUDA_VISIBLE_DEVICES="$gpu_id" protenix msa --input "$json_file" --out_dir "$MSA_DIR" > "$OUTPUT_DIR/${base_name}_msa.log" 2>&1; then
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
    CUDA_VISIBLE_DEVICES="$gpu_id" protenix pred \
        -i "$PREDICT_INPUT" \
        -o "$pred_dir" \
        -n "protenix_mini_default_v0.5.0" \
        --use_msa "$USE_MSA" \
        --use_default_params true \
        > "$OUTPUT_DIR/${base_name}_pred.log" 2>&1

    if [ $? -ne 0 ]; then
        log_ts "WARNING: protenix pred failed for $base_name. Check $OUTPUT_DIR/${base_name}_pred.log"
    else
        log_ts "Completed $base_name. Outputs saved to $pred_dir/"
    fi
}

if [ "$NUM_GPUS" -le 1 ]; then
    for json_file in "$JSON_DIR"/*.json; do
        [ -f "$json_file" ] || continue
        run_single "$json_file" "0"
    done
else
    log_ts "Multi-GPU mode: distributing across $NUM_GPUS GPUs"
    gpu_idx=0
    pids=()
    for json_file in "$JSON_DIR"/*.json; do
        [ -f "$json_file" ] || continue
        gpu_id=$((gpu_idx % NUM_GPUS))
        run_single "$json_file" "$gpu_id" &
        pids+=($!)
        gpu_idx=$((gpu_idx + 1))

        # Limit concurrent jobs to NUM_GPUS
        if [ ${#pids[@]} -ge "$NUM_GPUS" ]; then
            wait "${pids[0]}"
            pids=("${pids[@]:1}")
        fi
    done
    # Wait for remaining jobs
    for pid in "${pids[@]}"; do
        wait "$pid"
    done
fi

log_ts "Phase 1 High-Throughput Screening Complete."
