#!/bin/bash
# --- Option 4: Two-Phase MSA Workflow ---
# After Phase 1 screening (no MSA), re-run top N baselines WITH MSA for higher-quality structures.
# Baselines are ranked by Phase 1 Protenix ipTM score (from *_summary*.json); best first.
# Evolution then uses these MSA-enhanced structures as seeds.
#
# Usage: ./02b_rerun_top_with_msa.sh [N]
#   N = number of top baselines to re-run (default: 5, matches NUM_INITIAL_LINEAGES)
#
# Run AFTER 02_run_screening.sh completes. Expect ~30-50 min per baseline (remote MSA server).
# Ranks all baselines by ipTM. Validation filtering happens before evolution only.

log_ts() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

JSON_DIR="../jsons"
OUTPUT_DIR="../outputs/phase1_screening"
METADATA_FILE="../metadata/variant_domain_metadata.json"
TOP_N="${1:-5}"

if [ ! -f "$METADATA_FILE" ]; then
    log_ts "ERROR: Metadata not found at $METADATA_FILE. Run 01_parse_and_annotate.py first."
    exit 1
fi

# Discover *_summary*.json under Phase 1, load ipTM scores, sort by score (best first), take top N
log_ts "Ranking baselines by Phase 1 ipTM score..."
candidates=()
while IFS= read -r bid; do
    [ -z "$bid" ] && continue
    candidates+=("$bid")
done < <(python3 -c "
import json
import os
import sys
import glob

output_dir = os.path.abspath('$OUTPUT_DIR')
metadata_file = '$METADATA_FILE'
top_n = int('$TOP_N')

with open(metadata_file) as f:
    meta = json.load(f)

summaries = glob.glob(os.path.join(output_dir, '**', '*_summary*.json'), recursive=True)
baseline_scores = {}

for path in summaries:
    try:
        rel = os.path.relpath(path, output_dir)
    except ValueError:
        continue
    parts = rel.replace(chr(92), '/').split('/')
    if not parts or not parts[0].endswith('_pred'):
        continue
    bid = parts[0][:-5]
    if bid not in meta:
        continue
    try:
        with open(path) as f:
            data = json.load(f)
        score = float(data.get('iptm', data.get('ranking_score', 0.0)))
        if bid not in baseline_scores or score > baseline_scores[bid]:
            baseline_scores[bid] = score
    except (json.JSONDecodeError, KeyError, ValueError):
        continue

sorted_bids = sorted(baseline_scores.keys(), key=lambda b: baseline_scores[b], reverse=True)
for i, bid in enumerate(sorted_bids[:top_n], 1):
    print(f'  #{i} {bid} (ipTM={baseline_scores[bid]:.3f})', file=sys.stderr)
    print(bid)
")

if [ ${#candidates[@]} -eq 0 ]; then
    log_ts "ERROR: No Phase 1 summary JSONs found. Run 02_run_screening.sh first."
    exit 1
fi

count=${#candidates[@]}
log_ts "Re-running top $count baseline(s) by ipTM score with MSA (~30-50 min each)..."

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
