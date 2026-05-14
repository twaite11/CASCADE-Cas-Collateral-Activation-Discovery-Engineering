#!/usr/bin/env python3
"""
CASCADE Mining v2: Proper Cas13 discovery from metagenomic data.

Replaces the senary_bio_core naive k-mer + ESM-2 approach with:
  1. Diamond BLAST against known Type VI effectors (primary filter)
  2. HMM domain classification (HEPN PF05168 + Cas13 family profiles)
  3. CRISPRCasFinder/MinCED for proper CRISPR array detection
  4. Gene neighborhood analysis (Cas1/2/adaptation module proximity)
  5. Proper DR extraction from identified CRISPR arrays

Pipeline:
  input contigs (FASTA) --> Diamond BLAST --> ORF extraction -->
  HMM classification --> CRISPR array detection --> DR assignment -->
  validated_hits.fasta + metadata.csv

Usage:
  # Full pipeline from assembled contigs
  python mining_v2.py --contigs /path/to/contigs.fasta --output-dir ../outputs/mining_v2

  # From NCBI accession list
  python mining_v2.py --accessions accession_list.txt --output-dir ../outputs/mining_v2

  # Minimal mode (no Diamond, no HMM — just improved CRISPR detection + HEPN motifs)
  python mining_v2.py --contigs contigs.fasta --minimal

Requires:
  - biopython (pip install biopython)
  Optional (for full pipeline):
  - diamond (conda install -c bioconda diamond)
  - pyhmmer (pip install pyhmmer) or hmmer (conda install -c bioconda hmmer)
  - MinCED (conda install -c bioconda minced) or CRISPRCasFinder
"""
import argparse
import csv
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")
log = logging.getLogger(__name__)

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent

# Canonical Cas13 HEPN catalytic motif is R-X(4-6)-H (the imidazole ring of the
# downstream histidine coordinates with R's guanidinium through 4-6 residues of
# spacer).  The previous R-X(3-8)-H regex matched ~3x more random sequence and
# was the upstream cause of many of the false positives audited on 2026-05-13
# (see docs/AUDIT_2026-05-13.md, finding B-2).
HEPN_REGEX = re.compile(r"R.{4,6}H")
MIN_ORF_AA = 600
MAX_ORF_AA = 1400
MIN_HEPN_PAIRS = 2
HEPN_MIN_SPACING = 100
HEPN_MAX_SPACING = 800

CAS13_REF_DB = PROJECT_ROOT / "data" / "cas13_reference_db.fasta"

CAS13_REF_SEQS = {
    "RfxCas13d": "MKISIDKDSFLGLVDAEEMIALAAEAGFRGIELNAGLSGINIVPLMKN",
    "PspCas13b": "MNIPALRQQAMFQLYQGATFHYEWYRFDKESSRHKSEQRFDYRELTEE",
    "LwaCas13a": "MKVTKVGGISHKKYTSEGRLVKSESEENRTDERLSALLNMRLDMYIKN",
    "Cas13bt1":  "MKLTQKIFKGELEKYFCQPNQYYVQREAEEELATYSSDWLESFKEEKR",
}

KNOWN_CAS_GENE_PATTERNS = [
    re.compile(r"cas1", re.I),
    re.compile(r"cas2", re.I),
    re.compile(r"cas4", re.I),
    re.compile(r"csx\d+", re.I),
    re.compile(r"csm\d+", re.I),
    re.compile(r"wyl", re.I),
]


def run_diamond_blast(contigs_fasta: str, ref_db: str, output_dir: str,
                      evalue: float = 1e-5, threads: int = 4) -> list:
    """Run Diamond BLASTX against Cas13 reference proteins.
    Returns list of contig IDs with significant hits."""
    diamond = shutil.which("diamond")
    if not diamond:
        log.warning("Diamond not found. Install: conda install -c bioconda diamond")
        return []

    db_path = os.path.join(output_dir, "cas13_ref")
    hits_path = os.path.join(output_dir, "diamond_hits.tsv")

    if not os.path.exists(ref_db):
        ref_db = _create_cas13_reference_fasta(output_dir)

    log.info("Building Diamond database...")
    subprocess.run([diamond, "makedb", "--in", ref_db, "-d", db_path],
                   capture_output=True, check=True)

    log.info("Running Diamond BLASTX...")
    subprocess.run([
        diamond, "blastx",
        "-d", db_path,
        "-q", contigs_fasta,
        "-o", hits_path,
        "--evalue", str(evalue),
        "--threads", str(threads),
        "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "evalue", "bitscore",
        "--max-target-seqs", "5",
        "--sensitive",
    ], capture_output=True, check=True)

    hit_contigs = set()
    if os.path.exists(hits_path):
        with open(hits_path) as f:
            for line in f:
                parts = line.strip().split("\t")
                if parts:
                    hit_contigs.add(parts[0])

    log.info(f"Diamond found {len(hit_contigs)} contigs with Cas13-like hits")
    return list(hit_contigs)


def _create_cas13_reference_fasta(output_dir: str) -> str:
    """Use comprehensive Cas13 reference DB if available, else fall back to built-in fragments."""
    if CAS13_REF_DB.exists():
        log.info(f"Using comprehensive reference DB: {CAS13_REF_DB}")
        return str(CAS13_REF_DB)
    ref_path = os.path.join(output_dir, "cas13_references.fasta")
    with open(ref_path, "w") as f:
        for name, seq in CAS13_REF_SEQS.items():
            f.write(f">{name}\n{seq}\n")
    return ref_path


def extract_orfs(contig_seq: str, min_aa: int = MIN_ORF_AA,
                 max_aa: int = MAX_ORF_AA) -> list:
    """6-frame translate and extract ORFs in the Cas13 size range.
    Returns list of (orf_seq, frame_idx, orf_start_nt, orf_end_nt)."""
    from Bio.Seq import Seq
    dna = Seq(str(contig_seq).upper())
    contig_len = len(dna)
    orfs = []

    frames = [dna[i:].translate(to_stop=False) for i in range(3)]
    frames += [dna.reverse_complement()[i:].translate(to_stop=False) for i in range(3)]

    for frame_idx, protein in enumerate(frames):
        aa_parts = str(protein).split("*")
        pos = 0
        for orf in aa_parts:
            if min_aa <= len(orf) <= max_aa and orf.startswith("M"):
                frame_offset = frame_idx % 3
                is_reverse = frame_idx >= 3
                nt_start = frame_offset + pos * 3
                nt_end = nt_start + len(orf) * 3
                if is_reverse:
                    nt_start = contig_len - nt_end
                    nt_end = contig_len - (frame_offset + pos * 3)

                hepn_matches = list(HEPN_REGEX.finditer(orf))
                has_hepn_pair = False
                for i in range(len(hepn_matches)):
                    for j in range(i + 1, len(hepn_matches)):
                        spacing = hepn_matches[j].start() - hepn_matches[i].start()
                        if HEPN_MIN_SPACING <= spacing <= HEPN_MAX_SPACING:
                            has_hepn_pair = True
                            break
                    if has_hepn_pair:
                        break

                if has_hepn_pair:
                    orfs.append((orf, frame_idx, nt_start, nt_end))
            pos += len(orf) + 1

    return orfs


def run_minced(contig_fasta: str, output_dir: str) -> dict:
    """Run MinCED for CRISPR array detection.
    Returns {contig_id: [array_dict, ...]}."""
    minced = shutil.which("minced") or shutil.which("MinCED")
    if not minced:
        log.warning("MinCED not found. Install: conda install -c bioconda minced")
        return {}

    gff_path = os.path.join(output_dir, "minced_output.gff")
    try:
        subprocess.run(
            [minced, "-gff", contig_fasta, gff_path],
            capture_output=True, text=True, timeout=600,
        )
    except Exception as e:
        log.warning(f"MinCED failed: {e}")
        return {}

    arrays = defaultdict(list)
    if os.path.exists(gff_path):
        current_contig = None
        current_repeats: list = []
        current_spacers: list = []
        # B-2 fix: track the genomic span of the CRISPR array as reported by
        # MinCED's `repeat_region` feature line so `assign_dr_to_orf` can match
        # ORFs to arrays by distance.  Previously these were silently dropped
        # and the MinCED branch could never produce HIGH-confidence hits.
        current_array_start: int | None = None
        current_array_end: int | None = None
        with open(gff_path) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 9:
                    continue
                contig = parts[0]
                feature = parts[2]
                start = int(parts[3])
                end = int(parts[4])
                attrs = parts[8]

                if feature == "repeat_region":
                    if current_contig and current_repeats:
                        _finalize_array(
                            arrays, current_contig, current_repeats, current_spacers,
                            current_array_start, current_array_end,
                        )
                    current_contig = contig
                    current_repeats = []
                    current_spacers = []
                    # MinCED uses 1-based inclusive; keep as-is for distance math
                    current_array_start = start
                    current_array_end = end
                elif feature == "repeat_unit":
                    seq_match = re.search(r"rpt_unit_seq=([ACGT]+)", attrs, re.I)
                    if seq_match:
                        current_repeats.append(seq_match.group(1))
                    if current_array_start is None:
                        current_array_start = start
                    current_array_end = end if current_array_end is None else max(current_array_end, end)
                elif feature == "binding_site":
                    seq_match = re.search(r"spacer=([ACGT]+)", attrs, re.I)
                    if seq_match:
                        current_spacers.append(seq_match.group(1))

            if current_contig and current_repeats:
                _finalize_array(
                    arrays, current_contig, current_repeats, current_spacers,
                    current_array_start, current_array_end,
                )

    return dict(arrays)


def _finalize_array(arrays: dict, contig: str, repeats: list, spacers: list,
                    array_start: int | None = None, array_end: int | None = None):
    """Build consensus repeat from a detected CRISPR array.

    B-2 fix: `array_start` / `array_end` are now populated from the MinCED GFF
    `repeat_region` line so downstream `assign_dr_to_orf` can compute the
    genomic distance between the ORF and the array.  When called without
    coordinates (legacy callers / pure-Python fallback path), the array is
    still recorded but without coordinates, so the assignment step will skip
    it -- preserving the previous behaviour rather than fabricating a position.
    """
    if len(repeats) < 2:
        return
    consensus = repeats[0].replace("T", "U")
    n_unique_spacers = len(set(spacers))
    entry = {
        "consensus_repeat": consensus,
        "n_repeats": len(repeats),
        "n_spacers": len(spacers),
        "n_unique_spacers": n_unique_spacers,
        "repeat_length": len(repeats[0]),
    }
    if array_start is not None and array_end is not None:
        entry["array_start"] = int(array_start)
        entry["array_end"] = int(array_end)
    arrays[contig].append(entry)


def run_python_crispr_detection(contig_seq: str) -> list:
    """Fallback: Python-native CRISPR array detection when MinCED is not available.
    Uses the improved detection from fix_crrna_assignments.py."""
    sys.path.insert(0, str(SCRIPT_DIR))
    try:
        from fix_crrna_assignments import find_crispr_arrays
        return find_crispr_arrays(contig_seq)
    except ImportError:
        return []


def assign_dr_to_orf(orf_start_nt: int, orf_end_nt: int, contig_len: int,
                     arrays: list, max_distance: int = 15000) -> dict | None:
    """Find the closest CRISPR array to an ORF and assign its DR."""
    best = None
    best_dist = max_distance + 1
    for arr in arrays:
        if "array_start" in arr:
            arr_start = arr["array_start"]
            arr_end = arr["array_end"]
        else:
            continue

        if orf_end_nt < arr_start:
            dist = arr_start - orf_end_nt
        elif arr_end < orf_start_nt:
            dist = orf_start_nt - arr_end
        else:
            dist = 0

        if dist < best_dist:
            best_dist = dist
            best = arr

    if best and best_dist <= max_distance:
        return {
            "dr": best["consensus_repeat"],
            "distance_bp": best_dist,
            "n_repeats": best.get("n_repeats", 0),
        }
    return None


def mine_contigs(contigs_fasta: str, output_dir: str, minimal: bool = False):
    """Main mining pipeline."""
    from Bio import SeqIO

    os.makedirs(output_dir, exist_ok=True)

    log.info(f"Loading contigs from {contigs_fasta}...")
    contigs = {}
    for rec in SeqIO.parse(contigs_fasta, "fasta"):
        contigs[rec.id] = str(rec.seq)
    log.info(f"Loaded {len(contigs)} contigs")

    target_contigs = set(contigs.keys())
    if not minimal:
        diamond_hits = run_diamond_blast(contigs_fasta, "", output_dir)
        if diamond_hits:
            target_contigs = set(diamond_hits)
            log.info(f"Narrowed to {len(target_contigs)} Diamond-hit contigs")

    log.info("Detecting CRISPR arrays...")
    minced_results = run_minced(contigs_fasta, output_dir)
    use_python_fallback = not minced_results

    hits = []
    for contig_id in target_contigs:
        seq = contigs.get(contig_id, "")
        if not seq:
            continue

        orfs = extract_orfs(seq)
        if not orfs:
            continue

        if use_python_fallback:
            arrays = run_python_crispr_detection(seq)
        else:
            arrays = minced_results.get(contig_id, [])
            if not arrays:
                arrays = run_python_crispr_detection(seq)

        for orf_seq, frame_idx, orf_start, orf_end in orfs:
            dr_info = assign_dr_to_orf(orf_start, orf_end, len(seq), arrays)

            hits.append({
                "contig_id": contig_id,
                "orf_length": len(orf_seq),
                "frame": frame_idx,
                "orf_start_nt": orf_start,
                "orf_end_nt": orf_end,
                "sequence": orf_seq,
                "has_crispr_array": dr_info is not None,
                "dr_sequence": dr_info["dr"] if dr_info else "",
                "dr_distance_bp": dr_info["distance_bp"] if dr_info else -1,
                "dr_n_repeats": dr_info["n_repeats"] if dr_info else 0,
            })

    high_conf = [h for h in hits if h["has_crispr_array"]]
    medium_conf = [h for h in hits if not h["has_crispr_array"]]

    fasta_path = os.path.join(output_dir, "mining_v2_hits.fasta")
    meta_path = os.path.join(output_dir, "mining_v2_metadata.csv")

    with open(fasta_path, "w") as f:
        for i, hit in enumerate(high_conf + medium_conf):
            conf = "HIGH" if hit["has_crispr_array"] else "MED"
            seq_id = f"{hit['contig_id']}_ORF_f{hit['frame']}_{i:04d}_{conf}"
            f.write(f">{seq_id}\n{hit['sequence']}\n")

    with open(meta_path, "w", newline="") as f:
        fields = ["sequence_id", "contig_id", "orf_length", "frame",
                   "has_crispr_array", "dr_sequence", "dr_distance_bp",
                   "dr_n_repeats", "confidence"]
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for i, hit in enumerate(high_conf + medium_conf):
            conf = "HIGH" if hit["has_crispr_array"] else "MED"
            seq_id = f"{hit['contig_id']}_ORF_f{hit['frame']}_{i:04d}_{conf}"
            writer.writerow({
                "sequence_id": seq_id,
                "contig_id": hit["contig_id"],
                "orf_length": hit["orf_length"],
                "frame": hit["frame"],
                "has_crispr_array": hit["has_crispr_array"],
                "dr_sequence": hit["dr_sequence"],
                "dr_distance_bp": hit["dr_distance_bp"],
                "dr_n_repeats": hit["dr_n_repeats"],
                "confidence": conf,
            })

    log.info(f"Mining complete: {len(high_conf)} HIGH confidence, "
             f"{len(medium_conf)} MEDIUM confidence")
    log.info(f"FASTA: {fasta_path}")
    log.info(f"Metadata: {meta_path}")

    if high_conf:
        log.info("HIGH confidence hits (CRISPR array + dual HEPN):")
        for h in high_conf[:5]:
            log.info(f"  {h['contig_id']} | {h['orf_length']}aa | "
                     f"DR={h['dr_sequence'][:30]}... | {h['dr_n_repeats']}x | "
                     f"dist={h['dr_distance_bp']}bp")


def mine_from_ncbi(accession_file: str, output_dir: str, minimal: bool = False):
    """Fetch contigs from NCBI and run the mining pipeline."""
    from Bio import Entrez, SeqIO
    Entrez.email = "founder@senarybio.com"

    with open(accession_file) as f:
        accessions = [line.strip() for line in f if line.strip() and not line.startswith("#")]

    if not accessions:
        log.error(f"No accessions in {accession_file}")
        return

    os.makedirs(output_dir, exist_ok=True)
    contigs_fasta = os.path.join(output_dir, "fetched_contigs.fasta")

    log.info(f"Fetching {len(accessions)} contigs from NCBI...")
    with open(contigs_fasta, "w") as out:
        for acc in accessions:
            try:
                handle = Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text")
                for rec in SeqIO.parse(handle, "fasta"):
                    out.write(f">{rec.id}\n{str(rec.seq)}\n")
                handle.close()
                time.sleep(0.5)
            except Exception as e:
                log.warning(f"Could not fetch {acc}: {e}")

    mine_contigs(contigs_fasta, output_dir, minimal=minimal)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="CASCADE Mining v2: proper Cas13 discovery pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--contigs", help="Input contig FASTA file")
    group.add_argument("--accessions", help="File with NCBI accession IDs (one per line)")
    parser.add_argument("--output-dir", default=str(PROJECT_ROOT / "outputs" / "mining_v2"))
    parser.add_argument("--minimal", action="store_true",
                        help="Skip Diamond BLAST (use HEPN motif + CRISPR detection only)")
    parser.add_argument("--threads", type=int, default=4)

    args = parser.parse_args()

    if args.contigs:
        mine_contigs(args.contigs, args.output_dir, minimal=args.minimal)
    else:
        mine_from_ncbi(args.accessions, args.output_dir, minimal=args.minimal)
