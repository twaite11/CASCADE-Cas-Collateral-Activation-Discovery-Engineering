#!/bin/bash
# --- Option 4: Two-Phase MSA Workflow ---
# After Phase 1 screening (no MSA), re-run top N baselines WITH MSA for higher-quality PDBs.
# Evolution then uses these MSA-enhanced structures as seeds.
#
# Usage: ./02b_rerun_top_with_msa.sh [N]
#   N = number of top baselines to re-run (default: 5, matches NUM_INITIAL_LINEAGES)
#
# Run AFTER 02_run_screening.sh completes. Expect ~30-50 min per baseline (remote MSA server).

log_ts() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

JSON_DIR="../jsons"
OUTPUT_DIR="../outputs/phase1_screening"
METADATA_FILE="../metadata/variant_domain_metadata.json"
TOP_N="${1:-5}"

if [ ! -f "$METADATA_FILE" ]; then
    log_ts "ERROR: Metadata not found at $METADATA_FILE. Run 01_parse_and_annotate.py first."
    exit 1
fi

# Collect baseline IDs that have Phase 1 PDBs (same order as evolution_orchestrator)
candidates=()
while IFS= read -r bid; do
    [ -z "$bid" ] && continue
    if [ -d "$OUTPUT_DIR/${bid}_pred" ] && [ -n "$(find "$OUTPUT_DIR/${bid}_pred" -type f \( -name '*.pdb' -o -name '*.cif' \) 2>/dev/null | head -1)" ]; then
        candidates+=("$bid")
    fi
done < <(python3 -c "
import json
with open('$METADATA_FILE') as f:
    meta = json.load(f)
for k in list(meta.keys()):
    print(k)
")

if [ ${#candidates[@]} -eq 0 ]; then
    log_ts "ERROR: No Phase 1 PDBs found. Run 02_run_screening.sh first."
    exit 1
fi

# Take top N
count=$(( TOP_N < ${#candidates[@]} ? TOP_N : ${#candidates[@]} ))
log_ts "Re-running top $count of ${#candidates[@]} baselines with MSA (~30-50 min each)..."

for (( i=0; i<count; i++ )); do
    bid="${candidates[$i]}"
    json_file="$JSON_DIR/${bid}.json"
    if [ ! -f "$json_file" ]; then
        log_ts "WARNING: JSON not found for $bid, skipping."
        continue
    fi

    echo "=================================================="
    echo "MSA re-run $((i+1))/$count: $bid"
    echo "=================================================="

    MSA_DIR="$OUTPUT_DIR/${bid}_msa"
    mkdir -p "$MSA_DIR"

    log_ts "[1/2] Running MSA search (protenix msa) for $bid..."
    if ! protenix msa --input "$json_file" --out_dir "$MSA_DIR" \
        > "$OUTPUT_DIR/${bid}_msa.log" 2>&1; then
        log_ts "MSA failed for $bid. Keeping original no-MSA PDB. Check $OUTPUT_DIR/${bid}_msa.log"
        continue
    fi

    MSA_JSON="$MSA_DIR/$(basename "$json_file")"
    if [ ! -f "$MSA_JSON" ]; then
        log_ts "MSA output not found for $bid. Keeping original PDB."
        continue
    fi

    log_ts "[2/2] Running Protenix-Mini with MSA for $bid (replacing Phase 1 PDB)..."
    protenix predict \
        --input "$MSA_JSON" \
        --out_dir "$OUTPUT_DIR/${bid}_pred" \
        --model_name "protenix_mini_default_v0.5.0" \
        --use_msa true \
        --use_default_params true \
        > "$OUTPUT_DIR/${bid}_pred_msa.log" 2>&1

    log_ts "Completed MSA re-run for $bid."
done

log_ts "Phase 1b (MSA re-run for top $count) complete. Proceed to Phase 2: python evolution_orchestrator.py"
