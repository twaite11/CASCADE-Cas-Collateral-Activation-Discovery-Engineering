#!/bin/bash
# --- Option 4: Two-Phase MSA Workflow ---
# After Phase 1 screening (no MSA), re-run top N baselines WITH MSA for higher-quality structures.
# Baselines are ranked by Phase 1 ipTM score (from *_summary*.json or *_confidence*.json); best first.
# Supports cattle-prod (Rust, preferred) and protenix (Python, fallback). Set EVAL_CMD to override.
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

# --- Eval engine detection ---
if [ -n "$EVAL_CMD" ]; then
    EVAL_BIN="$EVAL_CMD"
elif command -v cattle-prod &>/dev/null; then
    EVAL_BIN="cattle-prod"
elif command -v protenix &>/dev/null; then
    EVAL_BIN="protenix"
else
    log_ts "ERROR: Neither cattle-prod nor protenix found on PATH. Set EVAL_CMD."
    exit 1
fi

if [[ "$EVAL_BIN" == *"cattle-prod"* ]]; then
    MINI_MODEL="cattle_prod_mini_default_v0.5.0"
    ENGINE_NAME="cattle-prod"
else
    MINI_MODEL="protenix_mini_default_v0.5.0"
    ENGINE_NAME="protenix"
fi
log_ts "Eval engine: $ENGINE_NAME ($EVAL_BIN)"

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

summaries = glob.glob(os.path.join(output_dir, '**', '*_summary*.json'), recursive=True) or glob.glob(os.path.join(output_dir, '**', '*_confidence*.json'), recursive=True)
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

    log_ts "[1/2] Running MSA search ($ENGINE_NAME msa) for $bid..."
    if ! "$EVAL_BIN" msa --input "$json_file" --out_dir "$MSA_DIR" \
        > "$OUTPUT_DIR/${bid}_msa.log" 2>&1; then
        log_ts "MSA failed for $bid. Keeping original no-MSA PDB. Check $OUTPUT_DIR/${bid}_msa.log"
        continue
    fi

    MSA_JSON="$JSON_DIR/${bid}-update-msa.json"
    if [ ! -f "$MSA_JSON" ]; then
        MSA_JSON="$MSA_DIR/$(basename "$json_file")"
    fi
    if [ ! -f "$MSA_JSON" ]; then
        log_ts "MSA output not found for $bid. Keeping original PDB."
        continue
    fi
    log_ts "Found MSA-updated JSON: $MSA_JSON"

    log_ts "[2/2] Running $ENGINE_NAME mini with MSA for $bid (replacing Phase 1 PDB)..."
    "$EVAL_BIN" pred \
        -i "$MSA_JSON" \
        -o "$OUTPUT_DIR/${bid}_pred" \
        -n "$MINI_MODEL" \
        --use_msa true \
        --use_default_params true \
        > "$OUTPUT_DIR/${bid}_pred_msa.log" 2>&1

    log_ts "Completed MSA re-run for $bid."
done

log_ts "Phase 1b (MSA re-run for top $count) complete. Proceed to Phase 2: python evolution_orchestrator.py"
