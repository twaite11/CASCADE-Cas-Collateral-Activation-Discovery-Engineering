#!/usr/bin/env python3
"""
CRISPR repeat validation: smarter k-mer selection + RNAfold structure filter.
Run right before evolution_orchestrator: python validate_crispr_repeats.py

Reads CSV files from data/mined_hits/, outputs:
  - outputs/repeat_validation_report.csv  (full report)
  - outputs/validated_baseline_ids.txt    (IDs passing threshold, one per line)

Use validated_baseline_ids.txt to restrict evolution_orchestrator to only baselines
with plausible CRISPR repeats. (02 screening processes all; filtering is at evolution.)
"""
import csv
import os
import re
import sys
from pathlib import Path

# Paths (relative to script location)
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
DATA_DIR = PROJECT_ROOT / "data" / "mined_hits"
OUTPUT_DIR = PROJECT_ROOT / "outputs"
REPORT_CSV = OUTPUT_DIR / "repeat_validation_report.csv"
VALIDATED_IDS_FILE = OUTPUT_DIR / "validated_baseline_ids.txt"

# Cas13 repeat length preferences
PREFERRED_MIN, PREFERRED_MAX = 28, 36
FALLBACK_MIN, FALLBACK_MAX = 23, 55


def is_simple_tandem_repeat(seq: str, unit_max: int = 6, min_repeats: int = 3) -> bool:
    """True if seq is dominated by a short repeating unit (e.g. CCGAAG repeated)."""
    seq = seq.upper()
    for unit_len in range(2, unit_max + 1):
        for i in range(len(seq) - unit_len * min_repeats + 1):
            unit = seq[i : i + unit_len]
            repeats = 1
            j = i + unit_len
            while j + unit_len <= len(seq) and seq[j : j + unit_len] == unit:
                repeats += 1
                j += unit_len
            if repeats >= min_repeats and repeats * unit_len >= len(seq) * 0.6:
                return True
    return False


def is_low_complexity(seq: str, max_single_frac: float = 0.5) -> bool:
    """True if >max_single_frac of bases are a single nucleotide."""
    from collections import Counter
    seq = seq.upper()
    counts = Counter(seq)
    if not counts:
        return True
    most_common = counts.most_common(1)[0][1]
    return (most_common / len(seq)) > max_single_frac


def select_best_repeat(k_mers: list) -> tuple:
    """Pick best repeat from k-mers. Returns (chosen_seq, reason)."""
    candidates = []
    for k in k_mers:
        k = k.strip().replace("T", "U").replace("t", "u")
        if not k or len(k) < FALLBACK_MIN or len(k) > FALLBACK_MAX:
            continue
        if is_simple_tandem_repeat(k) or is_low_complexity(k):
            continue
        candidates.append(k)

    if not candidates:
        return None, "no_valid_candidates"

    preferred = [c for c in candidates if PREFERRED_MIN <= len(c) <= PREFERRED_MAX]
    if preferred:
        chosen = max(preferred, key=len)
        return chosen, f"preferred_len_{len(chosen)}nt"
    chosen = max(candidates, key=len)
    return chosen, f"fallback_longest_{len(chosen)}nt"


def run_rnafold(seq: str):
    """Run RNAfold, return (structure, mfe) or (None, None) on failure."""
    try:
        import RNA  # pip install ViennaRNA  or  conda install -c bioconda viennarna
        structure, mfe = RNA.fold(seq)
        return structure, mfe
    except ImportError:
        pass
    try:
        import subprocess
        proc = subprocess.run(
            ["RNAfold", "--noPS"],
            input=seq.encode(),
            capture_output=True,
            timeout=5,
        )
        if proc.returncode != 0:
            return None, None
        line = proc.stdout.decode().strip().split("\n")[-1]
        parts = line.split()
        if not parts:
            return None, None
        # Robust MFE parse: RNAfold outputs "structure ( -2.60)" or "( 0.00)"
        # Anchor at end so we don't match "(....)" dots in structure; require digit
        mfe_match = re.search(r"\s*\([\s]*(-?\d+\.?\d*)\)\s*$", line)
        if mfe_match:
            mfe = float(mfe_match.group(1))
            structure = line[: mfe_match.start()].strip()
            return structure, mfe
    except Exception:
        pass
    return None, None


def has_plausible_stemloop(
    structure: str, mfe: float,
    min_bp: int = 4,
    max_mfe_per_nt: float = -0.05
) -> bool:
    """True if structure has stem (paired bases) and reasonable MFE."""
    if not structure:
        return False
    pairs = structure.count("(") + structure.count(")")
    if pairs < min_bp:
        return False
    n = len(structure)
    if n == 0:
        return False
    mfe_per_nt = mfe / n
    return mfe_per_nt <= max_mfe_per_nt


def load_validated_ids(path: Path = None) -> set:
    """Load validated baseline IDs from file. Returns empty set if file missing."""
    p = path or VALIDATED_IDS_FILE
    if not p.exists():
        return set()
    ids = set()
    with open(p) as f:
        for line in f:
            bid = line.strip()
            if bid:
                ids.add(bid)
    return ids


def main():
    # Handle CSV files with very large fields (e.g. long repeat_domains)
    csv.field_size_limit(min(2**31 - 1, sys.maxsize))

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    csv_files = list(DATA_DIR.glob("*_metadata.csv"))
    if not csv_files:
        print(f"No *_metadata.csv in {DATA_DIR}", file=sys.stderr)
        sys.exit(1)

    rows = []
    for csv_path in sorted(csv_files):
        with open(csv_path, encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                seq_id = row.get("sequence_id", "")
                raw = row.get("repeat_domains", "")
                k_mers = [x.strip() for x in raw.split("|") if x.strip()]

                chosen, reason = select_best_repeat(k_mers)
                structure, mfe = run_rnafold(chosen) if chosen else (None, None)
                structure_ok = has_plausible_stemloop(structure or "", mfe or 0) if chosen else False

                rows.append({
                    "sequence_id": seq_id,
                    "original_first_kmer": (k_mers[0][:50] if k_mers else ""),
                    "chosen_repeat": chosen or "",
                    "chosen_length": len(chosen) if chosen else 0,
                    "selection_reason": reason,
                    "structure": structure or "",
                    "mfe_kcal_mol": f"{mfe:.2f}" if mfe is not None else "",
                    "structure_ok": structure_ok,
                })

    with open(REPORT_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    validated = [r["sequence_id"] for r in rows if r["structure_ok"]]
    with open(VALIDATED_IDS_FILE, "w", encoding="utf-8") as f:
        for bid in validated:
            f.write(f"{bid}\n")

    print(f"Wrote {len(rows)} rows to {REPORT_CSV}")
    print(f"Passed structure filter: {len(validated)}/{len(rows)}")
    if len(validated) == 0 and any(r.get("chosen_repeat") for r in rows):
        print("  NOTE: Install ViennaRNA for structure filtering: pip install ViennaRNA")
    print(f"Validated IDs written to {VALIDATED_IDS_FILE}")
    if validated:
        print(f"  First 5: {', '.join(validated[:5])}")


if __name__ == "__main__":
    main()
