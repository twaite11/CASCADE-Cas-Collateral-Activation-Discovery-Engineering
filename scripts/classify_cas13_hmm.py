#!/usr/bin/env python3
"""
Classify putative Cas13 hits using HMM profiles (HEPN domain + Cas13 family).

Replaces the loose R.{4,6}H regex with proper domain-level classification using:
  1. Pfam PF05168 (HEPN domain) — the catalytic RNase domain
  2. Built-in Cas13 consensus HMMs (derived from known Cas13 a/b/d/X/Y subtypes)

Usage:
  python classify_cas13_hmm.py [--fasta-dir ../data/mined_hits] [--threshold 1e-5]

Requires: pyhmmer (pip install pyhmmer) OR hmmer (conda install -c bioconda hmmer)
"""
import csv
import json
import logging
import os
import re
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")
log = logging.getLogger(__name__)

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
DATA_DIR = PROJECT_ROOT / "data" / "mined_hits"
OUTPUT_DIR = PROJECT_ROOT / "outputs"
HMM_DIR = PROJECT_ROOT / "data" / "hmm"
METADATA_JSON = PROJECT_ROOT / "metadata" / "variant_domain_metadata.json"
VALIDATED_IDS = OUTPUT_DIR / "validated_baseline_ids.txt"
CLASSIFICATION_REPORT = OUTPUT_DIR / "cas13_classification_report.csv"

HEPN_REGEX = re.compile(r"R.{3,8}H")

HEPN_CONSENSUS_SEQS = [
    "RNYQKRGH",   # Cas13a HEPN1 consensus
    "RQHQNRGH",   # Cas13a HEPN2 consensus
    "RIYNQRGH",   # Cas13b consensus
    "RAAAVRH",    # Cas13d HEPN1 (RfxCas13d)
    "RQNHRRH",    # Cas13d HEPN2
    "RKYKQRGH",   # Cas13X
    "RRYRQRGH",   # Cas13Y
]

CAS13_DIAGNOSTIC_PATTERNS = {
    "cas13a": [
        re.compile(r"R.{4}H.{20,80}[LIVMF]{2}.{5,15}[DE].{5,20}[RK]"),
        re.compile(r"[LIVMF].{2}[DE].{10,30}R.{4}H"),
    ],
    "cas13b": [
        re.compile(r"R.{3,6}H.{100,400}R.{3,6}H"),
    ],
    "cas13d": [
        re.compile(r"R[A-Z]{3,5}H.{150,350}R[A-Z]{3,5}H"),
    ],
}


def classify_by_hepn_motifs(protein_seq: str) -> dict:
    """Classify a protein by its HEPN motif pattern and spacing.
    Returns dict with classification results."""
    seq = str(protein_seq).upper()
    matches = list(HEPN_REGEX.finditer(seq))

    result = {
        "n_hepn_motifs": len(matches),
        "hepn_positions": [(m.start(), m.end(), m.group()) for m in matches],
        "has_dual_hepn": False,
        "hepn_spacing": None,
        "subtype_guess": "unknown",
        "confidence": "low",
        "motif_quality": 0.0,
    }

    if len(matches) < 2:
        return result

    best_pair = None
    best_score = -1
    for i in range(len(matches)):
        for j in range(i + 1, len(matches)):
            spacing = matches[j].start() - matches[i].start()
            if 100 <= spacing <= 800:
                motif1 = matches[i].group()
                motif2 = matches[j].group()
                score = _score_hepn_pair(motif1, motif2, spacing)
                if score > best_score:
                    best_score = score
                    best_pair = (matches[i], matches[j], spacing, score)

    if best_pair is None:
        return result

    m1, m2, spacing, score = best_pair
    result["has_dual_hepn"] = True
    result["hepn_spacing"] = spacing
    result["motif_quality"] = score

    if score >= 0.6:
        result["confidence"] = "high"
    elif score >= 0.3:
        result["confidence"] = "medium"

    for subtype, patterns in CAS13_DIAGNOSTIC_PATTERNS.items():
        for pat in patterns:
            if pat.search(seq):
                result["subtype_guess"] = subtype
                break

    return result


def _score_hepn_pair(motif1: str, motif2: str, spacing: int) -> float:
    """Score a HEPN motif pair for Cas13 likelihood (0-1).
    Dual HEPN with Cas13-typical spacing (150-500aa) is the primary signal.
    Consensus motif similarity is a bonus, not a gate."""
    score = 0.0

    if 200 <= spacing <= 400:
        score += 0.35
    elif 150 <= spacing <= 500:
        score += 0.25
    elif 100 <= spacing <= 800:
        score += 0.10

    for motif in [motif1, motif2]:
        if motif[0] != "R" or motif[-1] != "H":
            continue
        inner = motif[1:-1]
        if 3 <= len(inner) <= 7:
            score += 0.10
        best_sim = 0.0
        for consensus in HEPN_CONSENSUS_SEQS:
            cons_inner = consensus[1:-1]
            matches = sum(1 for a, b in zip(inner, cons_inner) if a == b)
            sim = matches / max(len(inner), len(cons_inner), 1)
            best_sim = max(best_sim, sim)
        score += best_sim * 0.15

    return min(score, 1.0)


def try_pyhmmer_scan(fasta_path: str, hmm_path: str, evalue: float = 1e-5) -> dict:
    """Use pyhmmer for HMM scanning. Returns {seq_id: [hit_info, ...]}."""
    try:
        import pyhmmer
    except ImportError:
        return None

    results = defaultdict(list)
    try:
        with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
            hmms = list(hmm_file)
        with pyhmmer.easel.SequenceFile(fasta_path, digital=True,
                                         alphabet=pyhmmer.easel.Alphabet.amino()) as seq_file:
            seqs = list(seq_file)

        for hits in pyhmmer.hmmsearch(hmms, seqs, E=evalue):
            for hit in hits:
                if hit.included:
                    results[hit.name.decode()].append({
                        "hmm": hits.query_name.decode(),
                        "evalue": hit.evalue,
                        "score": hit.score,
                    })
    except Exception as e:
        log.warning(f"pyhmmer scan failed: {e}")
        return None

    return dict(results)


def try_hmmsearch_subprocess(fasta_path: str, hmm_path: str, evalue: float = 1e-5) -> dict:
    """Fall back to hmmsearch CLI if pyhmmer not available."""
    import shutil
    hmmsearch = shutil.which("hmmsearch")
    if not hmmsearch:
        return None

    results = defaultdict(list)
    try:
        import subprocess
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tbl", delete=False) as tbl:
            tbl_path = tbl.name

        proc = subprocess.run(
            [hmmsearch, "--tblout", tbl_path, "-E", str(evalue), hmm_path, fasta_path],
            capture_output=True, text=True, timeout=300,
        )
        if proc.returncode == 0 and os.path.exists(tbl_path):
            with open(tbl_path) as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    parts = line.split()
                    if len(parts) >= 5:
                        seq_id = parts[0]
                        hmm_name = parts[2]
                        e_val = float(parts[4])
                        score_val = float(parts[5])
                        results[seq_id].append({
                            "hmm": hmm_name,
                            "evalue": e_val,
                            "score": score_val,
                        })
        os.unlink(tbl_path)
    except Exception as e:
        log.warning(f"hmmsearch failed: {e}")
        return None

    return dict(results)


def load_sequences_from_fastas() -> dict:
    """Load all protein sequences from mined_hits FASTA files."""
    sequences = {}
    for fasta_path in sorted(DATA_DIR.glob("*.fasta")):
        current_id = ""
        current_seq = []
        with open(fasta_path) as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_id and current_seq:
                        sequences[current_id] = "".join(current_seq)
                    current_id = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)
            if current_id and current_seq:
                sequences[current_id] = "".join(current_seq)
    return sequences


def classify_all(evalue_threshold: float = 1e-5):
    """Run classification on all mined hits."""
    sequences = load_sequences_from_fastas()
    if not sequences:
        log.error(f"No FASTA files found in {DATA_DIR}")
        return

    log.info(f"Loaded {len(sequences)} protein sequences")

    hmm_results = None
    hmm_files = list(HMM_DIR.glob("*.hmm")) if HMM_DIR.exists() else []
    if hmm_files:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tmp:
            for sid, seq in sequences.items():
                tmp.write(f">{sid}\n{seq}\n")
            tmp_fasta = tmp.name

        for hmm_path in hmm_files:
            log.info(f"Scanning with HMM: {hmm_path.name}")
            result = try_pyhmmer_scan(tmp_fasta, str(hmm_path), evalue_threshold)
            if result is None:
                result = try_hmmsearch_subprocess(tmp_fasta, str(hmm_path), evalue_threshold)
            if result:
                if hmm_results is None:
                    hmm_results = {}
                for sid, hits in result.items():
                    hmm_results.setdefault(sid, []).extend(hits)

        os.unlink(tmp_fasta)

    validated_ids = set()
    if VALIDATED_IDS.exists():
        with open(VALIDATED_IDS) as f:
            validated_ids = {line.strip() for line in f if line.strip()}

    report_rows = []
    passed = []
    for seq_id, seq in sequences.items():
        motif_result = classify_by_hepn_motifs(seq)

        hmm_hits = hmm_results.get(seq_id, []) if hmm_results else []
        hmm_status = "hit" if hmm_hits else ("no_hmm_db" if not hmm_files else "no_hit")
        best_hmm_evalue = min((h["evalue"] for h in hmm_hits), default=None)

        is_cas13 = (
            motif_result["has_dual_hepn"]
            and motif_result["motif_quality"] >= 0.15
            and motif_result["confidence"] in ("medium", "high")
        )

        if hmm_hits:
            is_cas13 = True

        report_rows.append({
            "sequence_id": seq_id,
            "n_hepn": motif_result["n_hepn_motifs"],
            "dual_hepn": motif_result["has_dual_hepn"],
            "spacing": motif_result["hepn_spacing"] or "",
            "motif_quality": f"{motif_result['motif_quality']:.3f}",
            "confidence": motif_result["confidence"],
            "subtype": motif_result["subtype_guess"],
            "hmm_status": hmm_status,
            "hmm_evalue": f"{best_hmm_evalue:.1e}" if best_hmm_evalue else "",
            "is_cas13": is_cas13,
            "in_validated_ids": seq_id in validated_ids,
        })

        if is_cas13:
            passed.append(seq_id)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(CLASSIFICATION_REPORT, "w", newline="", encoding="utf-8") as f:
        if report_rows:
            writer = csv.DictWriter(f, fieldnames=list(report_rows[0].keys()))
            writer.writeheader()
            writer.writerows(report_rows)

    if validated_ids:
        refined = validated_ids & set(passed)
        with open(VALIDATED_IDS, "w") as f:
            for bid in sorted(refined):
                f.write(f"{bid}\n")
        log.info(f"Refined validated IDs: {len(validated_ids)} -> {len(refined)} "
                 f"(removed {len(validated_ids) - len(refined)} non-Cas13)")

    total = len(report_rows)
    cas13_count = sum(1 for r in report_rows if r["is_cas13"])
    high_conf = sum(1 for r in report_rows if r["confidence"] == "high")
    log.info(f"Classification: {cas13_count}/{total} classified as Cas13 "
             f"({high_conf} high confidence)")
    log.info(f"Report: {CLASSIFICATION_REPORT}")

    _propagate_subtypes_to_metadata(report_rows)

    if not hmm_files:
        log.info("TIP: Download HEPN HMM profiles for stronger classification:")
        log.info(f"  mkdir -p {HMM_DIR}")
        log.info(f"  wget -O {HMM_DIR}/HEPN.hmm 'https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/PF05168?annotation=hmm'")


def _propagate_subtypes_to_metadata(report_rows: list):
    """Write subtype_guess into variant_domain_metadata.json so downstream
    tools (protenix_eval, evolution_orchestrator) can use subtype-aware crRNA
    assembly (orientation, spacer length)."""
    if not METADATA_JSON.exists():
        log.info(f"No metadata JSON at {METADATA_JSON}; skipping subtype propagation.")
        return

    with open(METADATA_JSON) as f:
        metadata = json.load(f)

    updated = 0
    for row in report_rows:
        sid = row["sequence_id"]
        subtype = row.get("subtype", "unknown")
        if sid in metadata:
            if metadata[sid].get("subtype") != subtype:
                metadata[sid]["subtype"] = subtype
                updated += 1

    if updated:
        with open(METADATA_JSON, "w") as f:
            json.dump(metadata, f, indent=2)
        log.info(f"Propagated subtype to {updated} entries in {METADATA_JSON}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Classify Cas13 hits using HMM profiles")
    parser.add_argument("--fasta-dir", default=str(DATA_DIR))
    parser.add_argument("--threshold", type=float, default=1e-5)
    args = parser.parse_args()
    DATA_DIR = Path(args.fasta_dir)
    classify_all(evalue_threshold=args.threshold)
