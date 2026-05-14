#!/usr/bin/env python3
"""
Fix crRNA assignments: replace naive k-mer "CRISPR repeats" (often tRNAs)
with properly detected CRISPR direct repeats.

Pipeline:
  1. Read mined_hits metadata CSVs to get NCBI accessions
  2. Fetch source contigs from NCBI Entrez
  3. Detect real CRISPR arrays using repeat-spacer structure analysis
  4. Filter out tRNA-like sequences (cloverleaf, too-stable structure)
  5. Re-assign DRs to each ORF based on proximity
  6. Update variant_domain_metadata.json and validated_baseline_ids.txt

Usage:
  python fix_crrna_assignments.py [--offline]

  --offline   Skip NCBI fetch; only re-validate existing repeat_domains
              using the improved filters (tRNA exclusion, array structure).

Requires: biopython (pip install biopython)
Optional: ViennaRNA (pip install ViennaRNA) for structure filtering
"""
import csv
import json
import logging
import os
import re
import sys
import time
from collections import defaultdict
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")
log = logging.getLogger(__name__)

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
DATA_DIR = PROJECT_ROOT / "data" / "mined_hits"
METADATA_DIR = PROJECT_ROOT / "metadata"
OUTPUT_DIR = PROJECT_ROOT / "outputs"
METADATA_JSON = METADATA_DIR / "variant_domain_metadata.json"
VALIDATED_IDS = OUTPUT_DIR / "validated_baseline_ids.txt"
REPORT_CSV = OUTPUT_DIR / "repeat_validation_report.csv"
CONTIGS_CACHE = OUTPUT_DIR / "contig_cache"

CRISPR_REPEAT_MIN = 23
CRISPR_REPEAT_MAX = 50
CRISPR_SPACER_MIN = 15
CRISPR_SPACER_MAX = 80
MIN_ARRAY_UNITS = 3

_DNA_COMP = str.maketrans("ACGT", "TGCA")
_RNA_COMP = str.maketrans("ACGU", "UGCA")


def _reverse_complement_rna(rna_seq: str) -> str:
    """Return the reverse complement of an RNA sequence."""
    return rna_seq.upper().translate(_RNA_COMP)[::-1]


def _pick_best_dr_orientation(rna_dr: str) -> str:
    """Given an RNA DR, return whichever orientation (forward or reverse
    complement) produces a better single-stem-loop typical of a CRISPR
    direct repeat handle.

    Cas13 DR handles fold into a single stem-loop with MFE typically in the
    -3 to -12 kcal/mol range.  The correct orientation consistently produces
    a cleaner, more stable stem-loop than the reverse complement."""
    rc = _reverse_complement_rna(rna_dr)

    try:
        import RNA
    except ImportError:
        return rna_dr

    def _stem_loop_score(seq: str) -> float:
        structure, mfe = RNA.fold(seq)
        if not structure or len(seq) == 0:
            return -999.0
        pairs = structure.count("(")
        mfe_per_nt = mfe / len(seq)
        in_range = -0.35 < mfe_per_nt < -0.05
        return pairs * (1.0 if in_range else 0.3) + (0 if in_range else -50)

    fwd_score = _stem_loop_score(rna_dr)
    rc_score = _stem_loop_score(rc)

    if rc_score > fwd_score:
        log.debug(f"DR orientation flipped: {rna_dr[:15]}… → RC (fwd={fwd_score:.1f}, rc={rc_score:.1f})")
        return rc
    return rna_dr


KNOWN_TRNA_SEEDS = [
    "GCGGGTGTAGCTCAG",
    "GGGCCCGTAGCTCAG",
    "GCGCCGCTGGTCTAG",
    "GGAGCGTAGTTCAAT",
    "GCCGAAATAAGCGG",
    "GCCCGGATAGCTCAG",
    "GGCCGGTTAGCTCAG",
    "GCGGGAGTGGCGAAA",
    "GGGGCTATAGCTCAG",
    "GCCGCCGTAGCTCAG",
]

csv.field_size_limit(min(2**31 - 1, sys.maxsize))


def is_trna_like(seq: str) -> bool:
    """Detect if a sequence looks like a tRNA gene fragment.
    tRNAs are ~73-93nt, have a cloverleaf structure, and contain
    conserved sequence motifs. CRISPR DRs are typically 23-50nt."""
    seq_upper = seq.upper().replace("U", "T")
    if len(seq_upper) < 20:
        return False
    for seed in KNOWN_TRNA_SEEDS:
        if seed in seq_upper or seq_upper in seed:
            return True
    try:
        gc = (seq_upper.count("G") + seq_upper.count("C")) / len(seq_upper)
        if gc > 0.70 and len(seq_upper) >= 28:
            cag_count = seq_upper.count("CAG") + seq_upper.count("CCA")
            if cag_count >= 2:
                return True
    except ZeroDivisionError:
        pass
    return False


_VIENNA_MISSING_WARNED = False


def has_crispr_structure(seq: str) -> bool:
    """Check if a DR candidate has plausible CRISPR repeat structure via RNAfold.
    CRISPR DRs typically have a single stem-loop, MFE between -2 and -15 kcal/mol.
    tRNAs have much stronger structure (MFE < -20 kcal/mol per 70nt).

    B-10 fix: when ViennaRNA isn't installed, the previous behaviour was to
    `return True` (= accept any sequence as a DR candidate), which silently
    disabled the entire structural filter and let tRNA fragments through.  We
    now fail-closed: return False once, warn once at module level, and treat
    every subsequent call the same way until ViennaRNA is installed.  Callers
    that want the legacy permissive behaviour must opt in explicitly by setting
    the CASCADE_VIENNARNA_FAIL_OPEN environment variable to a truthy value.
    """
    global _VIENNA_MISSING_WARNED
    try:
        import RNA
        structure, mfe = RNA.fold(seq)
    except ImportError:
        if not _VIENNA_MISSING_WARNED:
            _VIENNA_MISSING_WARNED = True
            log.warning(
                "ViennaRNA not installed -- has_crispr_structure() failing closed. "
                "Install `viennarna` (conda-forge) to enable structural filtering, "
                "or set CASCADE_VIENNARNA_FAIL_OPEN=1 to restore the legacy "
                "permissive behaviour (NOT recommended for production mining)."
            )
        if os.environ.get("CASCADE_VIENNARNA_FAIL_OPEN", "").lower() in {"1", "true", "yes"}:
            return True
        return False
    if not structure:
        return False
    n = len(seq)
    if n == 0:
        return False
    mfe_per_nt = mfe / n
    pairs = structure.count("(")
    if pairs < 3:
        return False
    if mfe_per_nt < -0.35:
        return False  # too stable — likely tRNA or rRNA
    if mfe_per_nt > -0.02:
        return False  # too weak — no secondary structure
    return True


def find_crispr_arrays(dna_seq: str, min_units: int = MIN_ARRAY_UNITS) -> list:
    """Detect CRISPR arrays in a DNA contig using repeat-spacer structure analysis.

    Returns list of dicts: {
        "consensus_repeat": str,
        "repeat_positions": [(start, end), ...],
        "spacers": [str, ...],
        "array_start": int,
        "array_end": int,
        "n_repeats": int,
    }
    """
    seq = str(dna_seq).upper()
    length = len(seq)
    if length < 200:
        return []

    arrays = []
    for repeat_len in range(CRISPR_REPEAT_MIN, min(CRISPR_REPEAT_MAX + 1, length // 2)):
        seen = defaultdict(list)
        for i in range(length - repeat_len + 1):
            kmer = seq[i : i + repeat_len]
            if "N" in kmer:
                continue
            seen[kmer].append(i)

        for kmer, positions in seen.items():
            if len(positions) < min_units:
                continue

            positions.sort()
            cluster = _find_regular_cluster(positions, repeat_len)
            if cluster is None or len(cluster) < min_units:
                continue

            spacers = []
            valid = True
            for idx in range(len(cluster) - 1):
                spacer_start = cluster[idx] + repeat_len
                spacer_end = cluster[idx + 1]
                spacer_len = spacer_end - spacer_start
                if spacer_len < CRISPR_SPACER_MIN or spacer_len > CRISPR_SPACER_MAX:
                    valid = False
                    break
                spacers.append(seq[spacer_start:spacer_end])

            if not valid or len(spacers) < min_units - 1:
                continue

            spacer_lengths = [len(s) for s in spacers]
            mean_sp = sum(spacer_lengths) / len(spacer_lengths)
            length_var = max(spacer_lengths) - min(spacer_lengths)
            if length_var > 10:
                continue

            unique_spacers = set(spacers)
            if len(unique_spacers) < len(spacers) * 0.5:
                continue

            rna_kmer = kmer.replace("T", "U")
            if is_trna_like(rna_kmer):
                continue

            best_dr = _pick_best_dr_orientation(rna_kmer)

            arrays.append({
                "consensus_repeat": best_dr,
                "repeat_positions": [(p, p + repeat_len) for p in cluster],
                "spacers": spacers,
                "array_start": cluster[0],
                "array_end": cluster[-1] + repeat_len,
                "n_repeats": len(cluster),
            })

    arrays = _deduplicate_arrays(arrays)
    return arrays


def _find_regular_cluster(positions: list, repeat_len: int) -> list | None:
    """Find the largest subset of positions forming a regular repeat-spacer array.
    Repeats should be spaced by repeat_len + spacer_len, where spacer_len is
    roughly consistent."""
    if len(positions) < 3:
        return None

    best_cluster = None
    for i in range(len(positions)):
        cluster = [positions[i]]
        for j in range(i + 1, len(positions)):
            gap = positions[j] - cluster[-1]
            expected_min = repeat_len + CRISPR_SPACER_MIN
            expected_max = repeat_len + CRISPR_SPACER_MAX
            if expected_min <= gap <= expected_max:
                cluster.append(positions[j])
        if len(cluster) >= 3:
            if best_cluster is None or len(cluster) > len(best_cluster):
                best_cluster = cluster

    return best_cluster


def _cas13_dr_length_bonus(dr_len: int) -> float:
    """Score bonus for DR lengths matching known Cas13 subtypes.
    Cas13a: 31-36nt, Cas13b: 36nt, Cas13d: ~30nt.
    Returns 0-2 bonus points."""
    if 28 <= dr_len <= 36:
        return 2.0
    if 23 <= dr_len <= 27:
        return 0.5
    return 0.0


def _array_sort_key(arr: dict) -> float:
    """Composite score for ranking CRISPR arrays: n_repeats + DR length bonus."""
    dr_len = len(arr.get("consensus_repeat", ""))
    return arr["n_repeats"] + _cas13_dr_length_bonus(dr_len)


def _deduplicate_arrays(arrays: list) -> list:
    """Remove overlapping arrays, keeping the one with highest composite score
    (repeat count + Cas13 DR length preference)."""
    if not arrays:
        return []
    arrays.sort(key=_array_sort_key, reverse=True)
    kept = []
    used_ranges = []
    for arr in arrays:
        overlap = False
        for (s, e) in used_ranges:
            if arr["array_start"] < e and arr["array_end"] > s:
                overlap = True
                break
        if not overlap:
            kept.append(arr)
            used_ranges.append((arr["array_start"], arr["array_end"]))
    return kept


def fetch_contig(ncbi_id: str, cache_dir: Path) -> str | None:
    """Fetch a contig from NCBI Entrez, caching locally."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_file = cache_dir / f"{ncbi_id}.fasta"
    if cache_file.exists():
        with open(cache_file) as f:
            lines = f.readlines()
            return "".join(l.strip() for l in lines if not l.startswith(">"))

    try:
        from Bio import Entrez, SeqIO
        Entrez.email = "founder@senarybio.com"
        handle = Entrez.efetch(db="nucleotide", id=ncbi_id, rettype="fasta", retmode="text")
        records = list(SeqIO.parse(handle, "fasta"))
        handle.close()
        if not records:
            return None
        seq = str(records[0].seq)
        with open(cache_file, "w") as f:
            f.write(f">{ncbi_id}\n{seq}\n")
        time.sleep(0.5)
        return seq
    except Exception as e:
        log.warning(f"Could not fetch {ncbi_id}: {e}")
        return None


def reassign_crrna_for_baselines(offline: bool = False):
    """Main pipeline: re-scan contigs, find real CRISPR arrays, reassign DRs."""
    csv_files = sorted(DATA_DIR.glob("*_metadata.csv"))
    if not csv_files:
        log.error(f"No metadata CSVs in {DATA_DIR}")
        return

    accession_to_seqids = defaultdict(list)
    seqid_to_old_repeats = {}
    for csv_path in csv_files:
        with open(csv_path, encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                seq_id = row.get("sequence_id", "")
                accession = row.get("sra_accession", "")
                repeats_raw = row.get("repeat_domains", "")
                if seq_id and accession:
                    accession_to_seqids[accession].append(seq_id)
                    seqid_to_old_repeats[seq_id] = repeats_raw

    log.info(f"Found {len(seqid_to_old_repeats)} ORFs from {len(accession_to_seqids)} source contigs")

    accession_arrays = {}
    if not offline:
        log.info("Fetching source contigs from NCBI and scanning for CRISPR arrays...")
        for acc_id in accession_to_seqids:
            contig_seq = fetch_contig(acc_id, CONTIGS_CACHE)
            if contig_seq:
                arrays = find_crispr_arrays(contig_seq)
                accession_arrays[acc_id] = arrays
                if arrays:
                    arr_desc = ", ".join(str(a["n_repeats"]) + "x" + str(len(a["consensus_repeat"])) + "nt" for a in arrays)
                    log.info(f"  {acc_id}: {len(arrays)} CRISPR array(s) found ({arr_desc})")
                else:
                    log.warning(f"  {acc_id}: NO genuine CRISPR arrays detected")
            else:
                log.warning(f"  {acc_id}: could not fetch contig")

    report_rows = []
    validated_ids = []

    for seq_id, old_repeats in seqid_to_old_repeats.items():
        accession = seq_id.rsplit("_ORF_", 1)[0] if "_ORF_" in seq_id else ""
        old_kmers = [k.strip() for k in old_repeats.split("|") if k.strip()]

        new_dr = None
        source = "none"
        rejection_reason = ""

        if accession in accession_arrays and accession_arrays[accession]:
            best_array = max(accession_arrays[accession], key=_array_sort_key)
            new_dr = best_array["consensus_repeat"]
            source = f"crispr_array_{best_array['n_repeats']}x"
        elif not offline and accession in accession_arrays:
            rejection_reason = "no_crispr_array_on_contig"
        else:
            new_dr = _salvage_from_old_kmers(old_kmers)
            if new_dr:
                source = "salvaged_from_kmers"
            else:
                rejection_reason = "no_valid_repeat"

        is_valid = new_dr is not None
        report_rows.append({
            "sequence_id": seq_id,
            "accession": accession,
            "old_first_kmer": (old_kmers[0][:50] if old_kmers else ""),
            "new_dr": new_dr or "",
            "new_dr_length": len(new_dr) if new_dr else 0,
            "source": source,
            "rejection_reason": rejection_reason,
            "is_valid": is_valid,
        })
        if is_valid:
            validated_ids.append(seq_id)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(REPORT_CSV, "w", newline="", encoding="utf-8") as f:
        if report_rows:
            writer = csv.DictWriter(f, fieldnames=list(report_rows[0].keys()))
            writer.writeheader()
            writer.writerows(report_rows)

    with open(VALIDATED_IDS, "w") as f:
        for bid in validated_ids:
            f.write(f"{bid}\n")

    _update_metadata_json(report_rows)

    total = len(report_rows)
    valid = len(validated_ids)
    from_array = sum(1 for r in report_rows if r["source"].startswith("crispr_array"))
    salvaged = sum(1 for r in report_rows if r["source"] == "salvaged_from_kmers")
    rejected = total - valid
    log.info(f"Results: {valid}/{total} validated "
             f"({from_array} from real arrays, {salvaged} salvaged, {rejected} rejected)")
    log.info(f"Report: {REPORT_CSV}")
    log.info(f"Validated IDs: {VALIDATED_IDS}")


def _salvage_from_old_kmers(kmers: list) -> str | None:
    """Try to find a valid CRISPR DR from the old sliding-window k-mers
    by applying tRNA and structure filters."""
    candidates = []
    seen = set()
    for k in kmers:
        k_rna = k.strip().replace("T", "U").replace("t", "u")
        if len(k_rna) < CRISPR_REPEAT_MIN or len(k_rna) > CRISPR_REPEAT_MAX:
            continue
        if k_rna in seen:
            continue
        seen.add(k_rna)
        if is_trna_like(k_rna):
            continue
        if not has_crispr_structure(k_rna):
            continue
        candidates.append(k_rna)

    if not candidates:
        return None
    preferred = [c for c in candidates if 28 <= len(c) <= 36]
    chosen = preferred[0] if preferred else candidates[0]
    return _pick_best_dr_orientation(chosen)


def _update_metadata_json(report_rows: list):
    """Update variant_domain_metadata.json with corrected crRNA repeats."""
    if not METADATA_JSON.exists():
        log.info(f"No metadata JSON at {METADATA_JSON} (generated on VPS). "
                 f"Skipping JSON update — run 01_parse_and_annotate.py to regenerate.")
        return

    with open(METADATA_JSON) as f:
        metadata = json.load(f)

    dr_lookup = {r["sequence_id"]: r["new_dr"] for r in report_rows if r["new_dr"]}
    updated = 0
    for seq_id, data in metadata.items():
        if seq_id in dr_lookup:
            old = data.get("crRNA_repeat_used", "")
            new = dr_lookup[seq_id]
            if old != new:
                data["crRNA_repeat_used"] = new
                updated += 1

    METADATA_DIR.mkdir(parents=True, exist_ok=True)
    with open(METADATA_JSON, "w") as f:
        json.dump(metadata, f, indent=2)
    log.info(f"Updated {updated} crRNA entries in {METADATA_JSON}")


if __name__ == "__main__":
    offline = "--offline" in sys.argv
    if offline:
        log.info("Running in OFFLINE mode (no NCBI fetch, re-validating existing k-mers)")
    reassign_crrna_for_baselines(offline=offline)
