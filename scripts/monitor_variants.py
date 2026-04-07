#!/usr/bin/env python3
"""
Live monitor: compare each variant's amino acid changes against the original baseline.
Run in a second terminal while the orchestrator is running.

Usage:
    cd /workspace/CASCADE/scripts
    python monitor_variants.py              # one-shot scan
    python monitor_variants.py --watch 30   # refresh every 30 seconds
"""
import argparse
import glob
import json
import os
import sys
import time
from collections import defaultdict

SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.join(SCRIPTS_DIR, "..")
JSONS_DIR = os.path.join(BASE_DIR, "jsons")
METADATA_FILE = os.path.join(BASE_DIR, "metadata", "variant_domain_metadata.json")
GEN_DIR = os.path.join(BASE_DIR, "outputs", "generation_queue")
RL_DATASET = os.path.join(BASE_DIR, "outputs", "rl_gym_data", "rl_training_dataset.jsonl")

_AA_PROPERTIES = {
    "G": "tiny",    "A": "small",   "V": "hydrophobic", "L": "hydrophobic",
    "I": "hydrophobic", "P": "rigid", "F": "aromatic", "W": "aromatic",
    "M": "hydrophobic", "S": "polar", "T": "polar",    "C": "special",
    "Y": "aromatic", "H": "positive", "D": "negative", "E": "negative",
    "N": "polar",    "Q": "polar",   "K": "positive",  "R": "positive",
}

_PROPERTY_GROUPS = {
    "hydrophobic": {"V", "L", "I", "M", "F", "W", "A"},
    "polar":       {"S", "T", "N", "Q", "Y", "C"},
    "positive":    {"K", "R", "H"},
    "negative":    {"D", "E"},
    "special":     {"G", "P"},
}


def classify_change(old_aa, new_aa):
    """Classify a substitution as conservative, semi-conservative, or radical."""
    if old_aa == new_aa:
        return "identical"
    old_prop = _AA_PROPERTIES.get(old_aa, "?")
    new_prop = _AA_PROPERTIES.get(new_aa, "?")
    if old_prop == new_prop:
        return "conservative"
    shared = False
    for group_aas in _PROPERTY_GROUPS.values():
        if old_aa in group_aas and new_aa in group_aas:
            shared = True
            break
    return "semi-conservative" if shared else "radical"


def read_fasta_seq(path):
    """Read first sequence from FASTA."""
    seq_lines = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if seq_lines:
                    break
                continue
            seq_lines.append(line.strip())
    return "".join(seq_lines)


def get_baseline_seq(baseline_id):
    """Load the original baseline sequence from the base JSON."""
    json_path = os.path.join(JSONS_DIR, f"{baseline_id}.json")
    if os.path.exists(json_path):
        with open(json_path) as f:
            data = json.load(f)
        ent = data[0]["sequences"][0]
        prot = ent.get("proteinChain", ent.get("protein", {}))
        return prot.get("sequence", "")
    return None


def get_domain_coords(baseline_id):
    """Load HEPN domain boundaries from metadata."""
    if not os.path.exists(METADATA_FILE):
        return None
    with open(METADATA_FILE) as f:
        meta = json.load(f)
    entry = meta.get(baseline_id)
    if not entry:
        return None
    h1 = entry["domains"]["HEPN1"]
    h2 = entry["domains"]["HEPN2"]
    return {
        "hepn1": (h1["start"], h1["end"]),
        "hepn2": (h2["start"], h2["end"]),
    }


def region_label(pos, coords):
    """Return which domain region a 0-based position falls in."""
    if coords is None:
        return "?"
    h1s, h1e = coords["hepn1"]
    h2s, h2e = coords["hepn2"]
    if pos < h1s:
        return "REC/N-term"
    elif h1s <= pos < h1e:
        return "HEPN1"
    elif h1e <= pos < h2s:
        return "linker"
    elif h2s <= pos < h2e:
        return "HEPN2"
    else:
        return "C-term"


def diff_sequences(baseline_seq, variant_seq, coords=None):
    """Compare two sequences and return list of mutation dicts."""
    mutations = []
    max_len = max(len(baseline_seq), len(variant_seq))
    b = baseline_seq.ljust(max_len, "-")
    v = variant_seq.ljust(max_len, "-")
    for i in range(max_len):
        if b[i] != v[i]:
            mutations.append({
                "pos": i + 1,
                "old": b[i],
                "new": v[i],
                "region": region_label(i, coords),
                "change_type": classify_change(b[i], v[i]),
            })
    return mutations


def find_variant_fastas():
    """Find all variant FASTAs across all worker/generation directories."""
    patterns = [
        os.path.join(GEN_DIR, "**", "L*.fasta"),
        os.path.join(GEN_DIR, "**", "*.fasta"),
    ]
    found = set()
    for pat in patterns:
        for p in glob.glob(pat, recursive=True):
            found.add(p)
    return sorted(found, key=os.path.getmtime)


def load_rl_dataset():
    """Load fitness records from the RL training dataset."""
    records = {}
    paths = [RL_DATASET]
    paths += glob.glob(os.path.join(BASE_DIR, "outputs", "rl_gym_data", "worker_*", "rl_training_dataset.jsonl"))
    for path in paths:
        if not os.path.exists(path):
            continue
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    rec = json.loads(line)
                    vid = rec.get("variant_id", "")
                    if vid:
                        records[vid] = rec
                except json.JSONDecodeError:
                    continue
    return records


def _build_lineage_map():
    """Build a mapping from lineage hash prefix (e.g. 'L08c120') to baseline ID.

    Scans PXDesign YAML filenames and generation directories to find the
    original baseline ID associated with each compact lineage tag.
    """
    tag_to_bid = {}
    # Method 1: PXDesign YAMLs are named {baseline_id}_pxdesign_input.yaml
    for yml in glob.glob(os.path.join(GEN_DIR, "**", "*_pxdesign_input.yaml"), recursive=True):
        basename = os.path.basename(yml).replace("_pxdesign_input.yaml", "")
        # The lineage tag is derived from the crRNA lookup ID (which is the baseline ID for initial runs)
        import hashlib
        digest = hashlib.sha1(basename.encode("utf-8")).hexdigest()[:6]
        tag = f"L{digest}"
        tag_to_bid[tag] = basename
    # Method 2: directories named after the baseline ID inside gen dirs
    for d in glob.glob(os.path.join(GEN_DIR, "**", "3174*"), recursive=True):
        if os.path.isdir(d):
            basename = os.path.basename(d)
            import hashlib
            digest = hashlib.sha1(basename.encode("utf-8")).hexdigest()[:6]
            tag_to_bid[f"L{digest}"] = basename
    return tag_to_bid


def infer_baseline_id(fasta_path, lineage_map=None):
    """Try to figure out which baseline a variant came from."""
    parts = fasta_path.replace("\\", "/").split("/")
    # Direct match: directory contains the baseline ID
    for p in parts:
        if "ORF_Score" in p:
            return p
    # Reverse-map from lineage tag in the filename
    name = os.path.basename(fasta_path).replace(".fasta", "")
    lineage_tag = name.split("_")[0]  # e.g. 'L08c120'
    if lineage_map and lineage_tag in lineage_map:
        return lineage_map[lineage_tag]
    return None


def print_report(baseline_cache, rl_data):
    """Print the full variant diff report."""
    fastas = find_variant_fastas()
    if not fastas:
        print("\n  No variant FASTAs found yet. Waiting for the orchestrator to generate designs...\n")
        return

    lineage_map = _build_lineage_map()

    # Group by generation
    by_gen = defaultdict(list)
    for fp in fastas:
        name = os.path.basename(fp).replace(".fasta", "")
        parts = name.split("_")
        gen_tag = next((p for p in parts if p.startswith("g")), "g??")
        by_gen[gen_tag].append(fp)

    print("\n" + "=" * 90)
    print(f"  CASCADE VARIANT MONITOR — {len(fastas)} variants across {len(by_gen)} generation(s)")
    print(f"  Scanned: {GEN_DIR}")
    print("=" * 90)

    region_mutation_counts = defaultdict(int)
    change_type_counts = defaultdict(int)
    total_mutations = 0
    total_variants = 0

    for gen_tag in sorted(by_gen.keys()):
        gen_fastas = by_gen[gen_tag]
        print(f"\n{'─' * 90}")
        print(f"  Generation {gen_tag} — {len(gen_fastas)} variant(s)")
        print(f"{'─' * 90}")

        for fp in gen_fastas:
            name = os.path.basename(fp).replace(".fasta", "")
            variant_seq = read_fasta_seq(fp)
            if not variant_seq:
                continue

            # Try to find the baseline
            baseline_id = infer_baseline_id(fp, lineage_map)
            if baseline_id and baseline_id not in baseline_cache:
                seq = get_baseline_seq(baseline_id)
                coords = get_domain_coords(baseline_id)
                if seq:
                    baseline_cache[baseline_id] = (seq, coords)

            # Also try RL dataset for lineage info
            rl_rec = rl_data.get(name, {})
            rl_baseline = rl_rec.get("baseline_id", "")
            if rl_baseline and rl_baseline not in baseline_cache:
                seq = get_baseline_seq(rl_baseline)
                coords = get_domain_coords(rl_baseline)
                if seq:
                    baseline_cache[rl_baseline] = (seq, coords)

            bid = baseline_id or rl_baseline
            if bid and bid in baseline_cache:
                baseline_seq, coords = baseline_cache[bid]
            else:
                print(f"\n  {name}:  [baseline not found — cannot diff]")
                continue

            mutations = diff_sequences(baseline_seq, variant_seq, coords)
            total_variants += 1
            total_mutations += len(mutations)

            # Fitness from RL dataset
            fitness_str = ""
            if name in rl_data:
                r = rl_data[name]
                fitness_str = (f"  fitness={r.get('fitness', '?'):.1f}"
                              f"  iptm={r.get('iptm_score', '?')}"
                              f"  off={r.get('off_distance', '?'):.1f}A"
                              f"  on={r.get('on_distance', '?'):.1f}A")

            if not mutations:
                print(f"\n  {name}:  IDENTICAL to baseline (0 mutations){fitness_str}")
                continue

            # Summarize by region
            region_summary = defaultdict(int)
            for m in mutations:
                region_summary[m["region"]] += 1
                region_mutation_counts[m["region"]] += 1
                change_type_counts[m["change_type"]] += 1

            region_str = ", ".join(f"{r}:{c}" for r, c in sorted(region_summary.items()))
            print(f"\n  {name}:  {len(mutations)} mutation(s) [{region_str}]{fitness_str}")

            # Print individual mutations (cap at 30 for readability)
            show = mutations[:30]
            for m in show:
                change_icon = {"conservative": "~", "semi-conservative": "≈", "radical": "!"}
                icon = change_icon.get(m["change_type"], "?")
                print(f"    {icon} pos {m['pos']:>4d}  {m['old']}→{m['new']}  "
                      f"({m['region']}, {m['change_type']})")
            if len(mutations) > 30:
                print(f"    ... and {len(mutations) - 30} more")

    # Summary
    print(f"\n{'=' * 90}")
    print(f"  SUMMARY: {total_variants} variants, {total_mutations} total mutations")
    if region_mutation_counts:
        print(f"  By region:      {dict(region_mutation_counts)}")
    if change_type_counts:
        print(f"  By change type: {dict(change_type_counts)}")
    if total_variants > 0:
        print(f"  Avg mutations/variant: {total_mutations / total_variants:.1f}")

    # Flag concerns
    hepn_muts = region_mutation_counts.get("HEPN1", 0) + region_mutation_counts.get("HEPN2", 0)
    linker_muts = region_mutation_counts.get("linker", 0)
    if hepn_muts > 0:
        print(f"\n  ⚠ WARNING: {hepn_muts} mutation(s) in catalytic HEPN domains — "
              "these should be frozen! Check fixed_positions config.")
    if total_mutations > 0 and linker_muts == 0:
        print(f"\n  ⚠ WARNING: No linker mutations — the RL bias may not be reaching linker positions.")
    if total_mutations == 0 and total_variants > 0:
        print(f"\n  ⚠ WARNING: All variants are identical to baseline — "
              "MPNN may not be designing new sequences (check HETATM/ATOM fix).")
    print("=" * 90 + "\n")


def main():
    parser = argparse.ArgumentParser(description="Monitor CASCADE variant mutations vs baseline")
    parser.add_argument("--watch", type=int, default=0,
                       help="Refresh interval in seconds (0 = one-shot)")
    args = parser.parse_args()

    baseline_cache = {}

    if args.watch > 0:
        print(f"Watching every {args.watch}s — Ctrl+C to stop\n")
        while True:
            try:
                os.system("clear" if os.name != "nt" else "cls")
                rl_data = load_rl_dataset()
                print_report(baseline_cache, rl_data)
                print(f"  [Refreshing in {args.watch}s... Ctrl+C to stop]")
                time.sleep(args.watch)
            except KeyboardInterrupt:
                print("\nStopped.")
                break
    else:
        rl_data = load_rl_dataset()
        print_report(baseline_cache, rl_data)


if __name__ == "__main__":
    main()
