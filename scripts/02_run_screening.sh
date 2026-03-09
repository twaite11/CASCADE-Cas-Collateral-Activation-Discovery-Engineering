#!/bin/bash

# --- RunPod Phase 1 Execution Script ---
# This script executes the high-throughput 'True Cas' test utilizing
# the lightweight Protenix-Mini model to drastically save on compute costs.
# Uses protenix msa (when available) for better quality; falls back to raw JSON.
#
# Now generates BOTH OFF-state and ON-state structures:
#   OFF (Enzyme + gRNA)        → baseline structure for PXDesign (hyper-stabilize dormant)
#   ON  (Enzyme + gRNA + RNA)  → evaluates conformational activation potential
#
# Processes all baselines. Validation filtering happens before evolution only.

log_ts() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

JSON_DIR="../jsons"
OUTPUT_DIR="../outputs/phase1_screening"
SKIP_MSA="${SKIP_MSA:-1}"

mkdir -p "$OUTPUT_DIR"

# Collect unique baseline IDs from JSON filenames: {id}_OFF.json and {id}_ON.json
declare -A BASELINES
for json_file in "$JSON_DIR"/*_OFF.json; do
    [ -f "$json_file" ] || continue
    base_name=$(basename "$json_file" _OFF.json)
    BASELINES["$base_name"]=1
done

if [ ${#BASELINES[@]} -eq 0 ]; then
    log_ts "No *_OFF.json files found in $JSON_DIR. Checking for legacy single-JSON format..."
    # Fallback: legacy format where each baseline has a single {id}.json (ternary)
    for json_file in "$JSON_DIR"/*.json; do
        [ -f "$json_file" ] || continue
        base_name=$(basename "$json_file" .json)
        # Skip if it's an _OFF or _ON file
        [[ "$base_name" == *_OFF ]] && continue
        [[ "$base_name" == *_ON ]] && continue
        BASELINES["$base_name"]=1
    done
fi

log_ts "Found ${#BASELINES[@]} baselines to screen."

for base_name in "${!BASELINES[@]}"; do
    echo "=================================================="
    echo "Processing Hit: $base_name"
    echo "=================================================="

    # --- OFF State: Enzyme + gRNA (dormant) ---
    OFF_JSON="$JSON_DIR/${base_name}_OFF.json"
    if [ -f "$OFF_JSON" ]; then
        log_ts "[OFF] Running Protenix-Mini for OFF state (Enzyme + gRNA)..."
        PREDICT_INPUT="$OFF_JSON"
        USE_MSA="false"

        if [ "$SKIP_MSA" != "1" ]; then
            MSA_DIR="$OUTPUT_DIR/${base_name}_OFF_msa"
            mkdir -p "$MSA_DIR"
            if protenix msa --input "$OFF_JSON" --out_dir "$MSA_DIR" > "$OUTPUT_DIR/${base_name}_OFF_msa.log" 2>&1; then
                MSA_JSON="$MSA_DIR/$(basename "$OFF_JSON")"
                if [ -f "$MSA_JSON" ]; then
                    PREDICT_INPUT="$MSA_JSON"
                    USE_MSA="true"
                fi
            fi
        fi

        protenix predict \
            --input "$PREDICT_INPUT" \
            --out_dir "$OUTPUT_DIR/${base_name}_OFF_pred" \
            --model_name "protenix_mini_default_v0.5.0" \
            --use_msa "$USE_MSA" \
            --use_default_params true \
            > "$OUTPUT_DIR/${base_name}_OFF_pred.log" 2>&1

        log_ts "[OFF] Completed $base_name OFF state."
    else
        log_ts "[OFF] No OFF JSON found for $base_name, skipping OFF-state screening."
    fi

    # --- ON State: Enzyme + gRNA + Target RNA (activated) ---
    ON_JSON="$JSON_DIR/${base_name}_ON.json"
    if [ -f "$ON_JSON" ]; then
        log_ts "[ON] Running Protenix-Mini for ON state (Enzyme + gRNA + Target RNA)..."
        PREDICT_INPUT="$ON_JSON"
        USE_MSA="false"

        if [ "$SKIP_MSA" != "1" ]; then
            MSA_DIR="$OUTPUT_DIR/${base_name}_ON_msa"
            mkdir -p "$MSA_DIR"
            if protenix msa --input "$ON_JSON" --out_dir "$MSA_DIR" > "$OUTPUT_DIR/${base_name}_ON_msa.log" 2>&1; then
                MSA_JSON="$MSA_DIR/$(basename "$ON_JSON")"
                if [ -f "$MSA_JSON" ]; then
                    PREDICT_INPUT="$MSA_JSON"
                    USE_MSA="true"
                fi
            fi
        fi

        protenix predict \
            --input "$PREDICT_INPUT" \
            --out_dir "$OUTPUT_DIR/${base_name}_ON_pred" \
            --model_name "protenix_mini_default_v0.5.0" \
            --use_msa "$USE_MSA" \
            --use_default_params true \
            > "$OUTPUT_DIR/${base_name}_ON_pred.log" 2>&1

        log_ts "[ON] Completed $base_name ON state."
    else
        # Fallback: legacy single-JSON format
        LEGACY_JSON="$JSON_DIR/${base_name}.json"
        if [ -f "$LEGACY_JSON" ]; then
            log_ts "[Legacy] Running Protenix-Mini on legacy ternary JSON for $base_name..."
            protenix predict \
                --input "$LEGACY_JSON" \
                --out_dir "$OUTPUT_DIR/${base_name}_pred" \
                --model_name "protenix_mini_default_v0.5.0" \
                --use_msa "false" \
                --use_default_params true \
                > "$OUTPUT_DIR/${base_name}_pred.log" 2>&1
            log_ts "[Legacy] Completed $base_name."
        fi
    fi

    log_ts "Completed screening for $base_name."
done

log_ts "Phase 1 High-Throughput Screening Complete."
