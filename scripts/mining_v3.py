#!/usr/bin/env python3
"""CASCADE Mining v3 -- strict, biologically-correct Cas13 discovery.

Differences vs mining_v2.py (which produced ~58% non-Cas13 hits per the
2026-05-13 audit, see docs/AUDIT_2026-05-13.md):

  1. Strict canonical HEPN motif:  R-X(4-6)-H  (was R-X(3-8)-H)
  2. Anti-signature filter:        rejects ORFs carrying diagnostic non-Cas13
                                   fold motifs (signal peptide, lipobox, GH3
                                   beta-glucosidase, methionine synthase MetH,
                                   metallo-beta-lactamase / RNase Z, AAA+
                                   Walker A ATPases, sugar isomerases, Cas9
                                   LGLDLGTNSIGW signature).
  3. Reciprocal validation:        Biopython PairwiseAligner local alignment
                                   against data/cas13_reference_db.fasta, with
                                   per-orientation k-mer prefilter; ORF must
                                   reach >= --min-identity over >= --min-cov
                                   block of the best reference.
  4. MinCED-coord DR assignment:   uses Python-native find_crispr_arrays from
                                   fix_crrna_assignments.py which reliably
                                   populates array_start/array_end (mining_v2
                                   MinCED path is bugged; see B-9).
  5. Per-ORF rejection_reason:     every dropped ORF is logged with why so we
                                   can audit recall.

External tools (diamond, MinCED, hmmer, prodigal) are OPTIONAL. The default
configuration runs end-to-end in pure Python and reproduces the audit's
diagnosis of `data/recoverable_inputs/campaign2/` (Bacteroides + Flavobacterium
+ Leptotrichia contigs).

Usage:
    python scripts/mining_v3.py --contigs path/to/contigs.fasta \
        --output-dir outputs/mining_v3
    python scripts/mining_v3.py --contigs ... --quick   # skip ref alignment

Requires:  biopython
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import re
import sys
import time
from collections import Counter
from dataclasses import asdict, dataclass, field
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
DEFAULT_REF_DB = PROJECT_ROOT / "data" / "cas13_reference_db.fasta"

# ---------------------------------------------------------------------------
# Strict HEPN motif (canonical Cas13 catalytic dyad)
# ---------------------------------------------------------------------------
HEPN_STRICT = re.compile(r"R.{4,6}H")
HEPN_MIN_SPACING = 100   # aa between the two HEPN motifs of a single Cas13
HEPN_MAX_SPACING = 900

MIN_ORF_AA = 600
MAX_ORF_AA = 1400

# ---------------------------------------------------------------------------
# Diagnostic non-Cas13 fold motifs (drives the anti-signature filter)
# ---------------------------------------------------------------------------
SIGNATURES: dict[str, re.Pattern] = {
    "SIGNAL_PEPTIDE":     re.compile(r"^M[A-Z]{0,5}[KRN]?[FILVMWAYC]{6,15}[ACSGT]"),
    "LIPOBOX":            re.compile(r"^M[A-Z]{4,25}[LVI][ASTG][AGS]C"),
    "GH3_KHFPGHGD":       re.compile(r"K[HF]FPGHGD"),
    "GH3_CATALYTIC":      re.compile(r"R[LMV]DS[IT][EH]LYP[YFR]"),
    "METH_NTERM":         re.compile(r"DGAMGTM"),
    "METH_CTERM":         re.compile(r"M[TS]GMDEVGRLF"),
    "MBL_RNASE_Z":        re.compile(r"L[ST]HLHSDHV"),
    "WALKER_A":           re.compile(r"G[GAS]G[GAS]GKT"),
    "WALKER_A_ALT":       re.compile(r"K[VG][TVS]GK[ST]"),
    "SUGAR_ISOMERASE":    re.compile(r"IGGT|GGTFVD|DEAGT"),
    "CAS9_RUVC":          re.compile(r"LGLDLGTNSIGW"),
}

# Hard-disqualifying signature names (any one of these = reject)
HARD_DISQUALIFY = {
    "GH3_KHFPGHGD", "GH3_CATALYTIC",
    "METH_NTERM", "METH_CTERM",         # either alone is suspicious; both = METH
    "MBL_RNASE_Z",
    "WALKER_A", "WALKER_A_ALT",
    "SUGAR_ISOMERASE",
    "CAS9_RUVC",
}

# Soft signatures: only disqualify if no canonical HEPN pair found
SOFT_DISQUALIFY = {"SIGNAL_PEPTIDE", "LIPOBOX"}


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------
@dataclass
class OrfHit:
    contig_id: str
    sequence: str
    frame: int
    orf_start_nt: int
    orf_end_nt: int
    hepn_pair_count: int = 0
    flags: list[str] = field(default_factory=list)
    best_ref_hit: str = ""
    best_ref_identity: float = 0.0
    best_ref_block_len: int = 0
    has_crispr_array: bool = False
    dr_sequence: str = ""
    dr_distance_bp: int = -1
    dr_n_repeats: int = 0
    rejection_reason: str = ""
    confidence: str = "DROPPED"


# ---------------------------------------------------------------------------
# ORF extraction (size-filtered 6-frame, M-start; same logic as v2 minus
# the buggy HEPN pre-filter that used the loose regex)
# ---------------------------------------------------------------------------
def extract_orfs(contig_seq: str,
                 min_aa: int = MIN_ORF_AA,
                 max_aa: int = MAX_ORF_AA) -> list[tuple[str, int, int, int]]:
    from Bio.Seq import Seq

    dna = Seq(str(contig_seq).upper())
    contig_len = len(dna)
    orfs: list[tuple[str, int, int, int]] = []
    frames = [dna[i:].translate(to_stop=False) for i in range(3)]
    rc = dna.reverse_complement()
    frames += [rc[i:].translate(to_stop=False) for i in range(3)]

    for frame_idx, protein in enumerate(frames):
        aa_parts = str(protein).split("*")
        pos = 0
        for orf in aa_parts:
            if min_aa <= len(orf) <= max_aa and orf.startswith("M"):
                offset = frame_idx % 3
                is_rev = frame_idx >= 3
                nt_start = offset + pos * 3
                nt_end = nt_start + len(orf) * 3
                if is_rev:
                    nt_start, nt_end = contig_len - nt_end, contig_len - nt_start
                orfs.append((orf, frame_idx, nt_start, nt_end))
            pos += len(orf) + 1
    return orfs


# ---------------------------------------------------------------------------
# HEPN motif detection (strict)
# ---------------------------------------------------------------------------
def hepn_pairs(seq: str,
               min_sep: int = HEPN_MIN_SPACING,
               max_sep: int = HEPN_MAX_SPACING) -> int:
    """Return the count of valid HEPN motif pairs (R-X(4-6)-H separated by
    min_sep..max_sep aa)."""
    hits = [m.start() for m in HEPN_STRICT.finditer(seq)]
    n_pairs = 0
    for i, a in enumerate(hits):
        for b in hits[i + 1:]:
            if min_sep <= (b - a) <= max_sep:
                n_pairs += 1
    return n_pairs


# ---------------------------------------------------------------------------
# Anti-signature filter
# ---------------------------------------------------------------------------
def detect_signatures(seq: str) -> list[str]:
    flags: list[str] = []
    for name, pat in SIGNATURES.items():
        if name in {"SIGNAL_PEPTIDE", "LIPOBOX"}:
            if pat.match(seq[:30]):
                flags.append(name)
        elif pat.search(seq):
            flags.append(name)
    return flags


def is_disqualified(flags: list[str], hepn_pair_count: int) -> str:
    """Return rejection reason or '' if the ORF passes."""
    hard_hits = [f for f in flags if f in HARD_DISQUALIFY]
    if "METH_NTERM" in flags and "METH_CTERM" in flags:
        return "METH_full_signature(N+C)"
    if hard_hits:
        return f"hard_signature:{','.join(hard_hits)}"
    soft_hits = [f for f in flags if f in SOFT_DISQUALIFY]
    if soft_hits and hepn_pair_count == 0:
        return f"soft_signature_no_hepn:{','.join(soft_hits)}"
    return ""


# ---------------------------------------------------------------------------
# Reciprocal validation (Biopython local PairwiseAligner)
# ---------------------------------------------------------------------------
def _load_aligner():
    from Bio.Align import PairwiseAligner, substitution_matrices

    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    return aligner


def _kmer_overlap(query: str, ref: str, k: int = 5) -> int:
    """Number of distinct k-mers shared between query and ref (cheap prefilter)."""
    qkmers = {query[i:i + k] for i in range(len(query) - k + 1)}
    rkmers = {ref[i:i + k] for i in range(len(ref) - k + 1)}
    return len(qkmers & rkmers)


def reciprocal_validate(orf: str,
                        ref_seqs: dict[str, str],
                        aligner,
                        min_identity: float = 0.30,
                        min_block_len: int = 80,
                        kmer_floor: int = 25) -> tuple[str, float, int]:
    """Local-align orf against each reference Cas13.

    Returns (best_ref_id, best_identity, best_block_len). All zero if nothing
    passes the kmer pre-filter (which means the ORF doesn't share even tiny
    fragments with any known Cas13 - a strong negative signal).
    """
    best = ("", 0.0, 0)
    for ref_id, ref_seq in ref_seqs.items():
        if _kmer_overlap(orf, ref_seq, k=5) < kmer_floor:
            continue
        try:
            alignments = aligner.align(orf, ref_seq)
            aln = alignments[0]
        except (ValueError, OverflowError):
            continue
        # Compute identity over the aligned block
        q_aligned, r_aligned = str(aln).splitlines()[0], str(aln).splitlines()[2]
        # Biopython renders alignment as 3 lines; safer to count via aln.counts()
        try:
            counts = aln.counts()
            matches = counts.identities
            block_len = counts.aligned
        except Exception:
            # Fallback: count from the rendered string
            mid = str(aln).splitlines()[1]
            matches = sum(1 for c in mid if c not in {" ", ".", "-", "X"})
            block_len = max(len(q_aligned.replace("-", "")), 1)
        identity = matches / block_len if block_len else 0.0
        if identity >= min_identity and block_len >= min_block_len:
            if identity * block_len > best[1] * best[2]:
                best = (ref_id, identity, block_len)
    return best


# ---------------------------------------------------------------------------
# CRISPR array detection (Python-native via fix_crrna_assignments)
# ---------------------------------------------------------------------------
def detect_crispr_arrays(contig_seq: str) -> list[dict]:
    sys.path.insert(0, str(SCRIPT_DIR))
    from fix_crrna_assignments import find_crispr_arrays
    return find_crispr_arrays(contig_seq)


def assign_dr_to_orf(orf_start: int, orf_end: int,
                     arrays: list[dict],
                     max_distance: int = 15000) -> dict | None:
    best = None
    best_dist = max_distance + 1
    for arr in arrays:
        a_start = arr.get("array_start")
        a_end = arr.get("array_end")
        if a_start is None or a_end is None:
            continue
        if orf_end < a_start:
            dist = a_start - orf_end
        elif a_end < orf_start:
            dist = orf_start - a_end
        else:
            dist = 0
        if dist < best_dist:
            best_dist = dist
            best = arr
    if best is not None and best_dist <= max_distance:
        return {
            "dr": best["consensus_repeat"],
            "distance_bp": best_dist,
            "n_repeats": best.get("n_repeats", 0),
        }
    return None


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------
def load_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    cur: str | None = None
    buf: list[str] = []
    with path.open() as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if cur is not None:
                    seqs[cur] = "".join(buf)
                cur = line[1:].split()[0]
                buf = []
            elif line:
                buf.append(line)
    if cur is not None:
        seqs[cur] = "".join(buf)
    return seqs


def mine(contigs_path: Path,
         output_dir: Path,
         ref_db_path: Path = DEFAULT_REF_DB,
         min_identity: float = 0.30,
         min_block_len: int = 80,
         skip_reciprocal: bool = False) -> list[OrfHit]:
    output_dir.mkdir(parents=True, exist_ok=True)

    contigs = load_fasta(contigs_path)
    log.info(f"Loaded {len(contigs)} contigs from {contigs_path}")

    ref_seqs = load_fasta(ref_db_path) if not skip_reciprocal else {}
    log.info(f"Loaded {len(ref_seqs)} reference Cas13 proteins from "
             f"{ref_db_path.name}" if ref_seqs else "Reciprocal validation: SKIPPED")

    aligner = _load_aligner() if ref_seqs else None
    rejection_counter: Counter[str] = Counter()
    all_hits: list[OrfHit] = []

    for contig_id, contig_seq in contigs.items():
        t0 = time.time()
        orfs = extract_orfs(contig_seq)
        if not orfs:
            log.info(f"{contig_id}: 0 ORFs in size range; skipping")
            continue

        arrays = detect_crispr_arrays(contig_seq)
        log.info(f"{contig_id}: {len(orfs)} ORFs in size range, "
                 f"{len(arrays)} CRISPR arrays")

        for orf_seq, frame, nt_s, nt_e in orfs:
            hit = OrfHit(
                contig_id=contig_id,
                sequence=orf_seq,
                frame=frame,
                orf_start_nt=nt_s,
                orf_end_nt=nt_e,
            )

            # 1. HEPN strict filter
            hit.hepn_pair_count = hepn_pairs(orf_seq)
            if hit.hepn_pair_count == 0:
                hit.rejection_reason = "no_canonical_R-X(4-6)-H_pair"
                rejection_counter[hit.rejection_reason] += 1
                all_hits.append(hit)
                continue

            # 2. Anti-signature filter
            hit.flags = detect_signatures(orf_seq)
            disqual = is_disqualified(hit.flags, hit.hepn_pair_count)
            if disqual:
                hit.rejection_reason = disqual
                rejection_counter[disqual] += 1
                all_hits.append(hit)
                continue

            # 3. Reciprocal validation (optional)
            if aligner is not None:
                ref_id, ident, block = reciprocal_validate(
                    orf_seq, ref_seqs, aligner,
                    min_identity=min_identity,
                    min_block_len=min_block_len,
                )
                hit.best_ref_hit = ref_id
                hit.best_ref_identity = round(ident, 3)
                hit.best_ref_block_len = block
                if not ref_id:
                    hit.rejection_reason = (
                        f"no_reciprocal_hit(>={min_identity:.0%}id over "
                        f"{min_block_len}aa)"
                    )
                    rejection_counter[hit.rejection_reason] += 1
                    all_hits.append(hit)
                    continue

            # 4. CRISPR array assignment
            dr = assign_dr_to_orf(nt_s, nt_e, arrays)
            if dr:
                hit.has_crispr_array = True
                hit.dr_sequence = dr["dr"]
                hit.dr_distance_bp = dr["distance_bp"]
                hit.dr_n_repeats = dr["n_repeats"]
                hit.confidence = "HIGH"
            else:
                hit.confidence = "MED"

            all_hits.append(hit)

        log.info(f"{contig_id}: processed in {time.time()-t0:.1f}s")

    # ---- Write outputs ----------------------------------------------------
    accepted = [h for h in all_hits if h.confidence != "DROPPED"]
    accepted.sort(key=lambda h: (h.confidence != "HIGH", -h.best_ref_identity))

    fasta_path = output_dir / "mining_v3_hits.fasta"
    meta_path = output_dir / "mining_v3_metadata.csv"
    rejected_path = output_dir / "mining_v3_rejected.csv"

    with fasta_path.open("w") as fh:
        for h in accepted:
            seq_id = (
                f"{h.contig_id}_ORF_f{h.frame}_{h.orf_start_nt:09d}_{h.confidence}"
            )
            fh.write(f">{seq_id}\n{h.sequence}\n")

    meta_fields = [
        "sequence_id", "contig_id", "frame", "orf_start_nt", "orf_end_nt",
        "orf_length", "hepn_pair_count", "flags",
        "best_ref_hit", "best_ref_identity", "best_ref_block_len",
        "has_crispr_array", "dr_sequence", "dr_distance_bp", "dr_n_repeats",
        "confidence",
    ]
    with meta_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=meta_fields)
        w.writeheader()
        for h in accepted:
            seq_id = (
                f"{h.contig_id}_ORF_f{h.frame}_{h.orf_start_nt:09d}_{h.confidence}"
            )
            w.writerow({
                "sequence_id": seq_id,
                "contig_id": h.contig_id,
                "frame": h.frame,
                "orf_start_nt": h.orf_start_nt,
                "orf_end_nt": h.orf_end_nt,
                "orf_length": len(h.sequence),
                "hepn_pair_count": h.hepn_pair_count,
                "flags": "|".join(h.flags),
                "best_ref_hit": h.best_ref_hit,
                "best_ref_identity": h.best_ref_identity,
                "best_ref_block_len": h.best_ref_block_len,
                "has_crispr_array": h.has_crispr_array,
                "dr_sequence": h.dr_sequence,
                "dr_distance_bp": h.dr_distance_bp,
                "dr_n_repeats": h.dr_n_repeats,
                "confidence": h.confidence,
            })

    rej_fields = [
        "contig_id", "frame", "orf_start_nt", "orf_length",
        "hepn_pair_count", "flags", "rejection_reason",
    ]
    with rejected_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=rej_fields)
        w.writeheader()
        for h in all_hits:
            if h.confidence != "DROPPED":
                continue
            w.writerow({
                "contig_id": h.contig_id,
                "frame": h.frame,
                "orf_start_nt": h.orf_start_nt,
                "orf_length": len(h.sequence),
                "hepn_pair_count": h.hepn_pair_count,
                "flags": "|".join(h.flags),
                "rejection_reason": h.rejection_reason,
            })

    summary = {
        "contigs_scanned": len(contigs),
        "orfs_evaluated": len(all_hits),
        "accepted": len(accepted),
        "high_confidence": sum(1 for h in accepted if h.confidence == "HIGH"),
        "medium_confidence": sum(1 for h in accepted if h.confidence == "MED"),
        "rejection_counts": dict(rejection_counter),
        "min_identity": min_identity,
        "min_block_len": min_block_len,
        "skip_reciprocal": skip_reciprocal,
    }
    (output_dir / "mining_v3_summary.json").write_text(
        json.dumps(summary, indent=2)
    )

    log.info("=" * 70)
    log.info(f"Mining v3 complete:")
    log.info(f"  contigs scanned:   {summary['contigs_scanned']}")
    log.info(f"  ORFs evaluated:    {summary['orfs_evaluated']}")
    log.info(f"  accepted:          {summary['accepted']} "
             f"(HIGH={summary['high_confidence']}, "
             f"MED={summary['medium_confidence']})")
    log.info(f"  rejection counts:")
    for reason, n in rejection_counter.most_common():
        log.info(f"    {n:5d}  {reason}")
    log.info(f"FASTA:    {fasta_path}")
    log.info(f"Metadata: {meta_path}")
    log.info(f"Rejected: {rejected_path}")
    return all_hits


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--contigs", required=True, type=Path,
                   help="Input contigs FASTA")
    p.add_argument("--output-dir", type=Path,
                   default=PROJECT_ROOT / "outputs" / "mining_v3")
    p.add_argument("--ref-db", type=Path, default=DEFAULT_REF_DB)
    p.add_argument("--min-identity", type=float, default=0.30,
                   help="Min identity over local block to accept reciprocal hit")
    p.add_argument("--min-block-len", type=int, default=80,
                   help="Min aligned block length (aa) for reciprocal validation")
    p.add_argument("--quick", action="store_true",
                   help="Skip reciprocal validation (motif + signature only)")
    args = p.parse_args()

    mine(
        contigs_path=args.contigs,
        output_dir=args.output_dir,
        ref_db_path=args.ref_db,
        min_identity=args.min_identity,
        min_block_len=args.min_block_len,
        skip_reciprocal=args.quick,
    )


if __name__ == "__main__":
    main()
