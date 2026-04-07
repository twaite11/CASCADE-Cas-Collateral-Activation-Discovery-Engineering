#!/usr/bin/env python3
"""
DIAMOND BLAST all candidate proteins against comprehensive Cas13 reference database.

Covers 8 characterized Cas13 proteins across 4 subtypes:
  - Cas13a: LshCas13a (1389aa), LseCas13a (1139aa), LbuCas13a (1159aa)
  - Cas13b: PbuCas13b (1127aa), BzCas13b (1116aa)
  - Cas13d: EsCas13d (954aa), RfxCas13d/CasRx (967aa)
  - Cas13bt: Cas13bt3 (775aa)

Usage:
  python blast_candidates.py                    # BLAST all mined hits
  python blast_candidates.py --evalue 0.01      # Relaxed e-value threshold
  python blast_candidates.py --mode blastx      # BLASTX DNA contigs (6-frame)
  python blast_candidates.py --contigs /tmp/all_contigs.fasta --mode blastx
"""
import argparse
import glob
import logging
import os
import shutil
import subprocess
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")
log = logging.getLogger(__name__)

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
DATA_DIR = PROJECT_ROOT / "data" / "mined_hits"
REF_DB = PROJECT_ROOT / "data" / "cas13_reference_db.fasta"
OUTPUT_DIR = PROJECT_ROOT / "outputs"


def combine_fastas(fasta_dir: Path, output_path: str) -> int:
    """Combine all FASTA files into one, return sequence count."""
    count = 0
    with open(output_path, "w") as out:
        for fp in sorted(fasta_dir.glob("*.fasta")):
            with open(fp) as f:
                for line in f:
                    out.write(line)
                    if line.startswith(">"):
                        count += 1
    return count


def run_blast(query_fasta: str, ref_fasta: str, output_tsv: str,
              mode: str = "blastp", evalue: float = 1e-3,
              threads: int = 4) -> list:
    diamond = shutil.which("diamond")
    if not diamond:
        log.error("DIAMOND not found. Install: pip install diamond-bio  OR  conda install -c bioconda diamond")
        sys.exit(1)

    db_path = "/tmp/cas13_comprehensive_db"
    log.info(f"Building DIAMOND database from {ref_fasta}...")
    subprocess.run([diamond, "makedb", "--in", ref_fasta, "-d", db_path],
                   capture_output=True, check=True)

    blast_cmd = mode if mode in ("blastp", "blastx") else "blastp"
    log.info(f"Running DIAMOND {blast_cmd} (evalue={evalue})...")
    cmd = [
        diamond, blast_cmd,
        "-d", db_path,
        "-q", query_fasta,
        "-o", output_tsv,
        "--outfmt", "6", "qseqid", "sseqid", "pident", "length",
        "evalue", "bitscore", "qlen", "slen",
        "--evalue", str(evalue),
        "--threads", str(threads),
        "--max-target-seqs", "5",
        "--very-sensitive",
    ]
    if mode == "blastx":
        cmd.extend(["--frameshift", "15"])

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log.error(f"DIAMOND failed: {result.stderr}")
        return []

    hits = []
    if os.path.exists(output_tsv):
        with open(output_tsv) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 8:
                    hits.append({
                        "query": parts[0],
                        "subject": parts[1],
                        "pident": float(parts[2]),
                        "length": int(parts[3]),
                        "evalue": float(parts[4]),
                        "bitscore": float(parts[5]),
                        "qlen": int(parts[6]),
                        "slen": int(parts[7]),
                    })
    return hits


def main():
    parser = argparse.ArgumentParser(description="BLAST candidates against comprehensive Cas13 DB")
    parser.add_argument("--evalue", type=float, default=1e-3)
    parser.add_argument("--mode", choices=["blastp", "blastx"], default="blastp",
                        help="blastp for protein queries, blastx for DNA contigs")
    parser.add_argument("--contigs", default=None,
                        help="Path to DNA contigs FASTA (for --mode blastx)")
    parser.add_argument("--query", default=None,
                        help="Path to protein FASTA (for --mode blastp)")
    parser.add_argument("--ref", default=str(REF_DB),
                        help="Path to reference FASTA")
    args = parser.parse_args()

    if not os.path.exists(args.ref):
        log.error(f"Reference DB not found: {args.ref}")
        sys.exit(1)

    if args.mode == "blastx":
        if not args.contigs:
            log.error("--contigs is required for blastx mode")
            sys.exit(1)
        query_fasta = args.contigs
        n_seqs = sum(1 for line in open(query_fasta) if line.startswith(">"))
    else:
        if args.query:
            query_fasta = args.query
            n_seqs = sum(1 for line in open(query_fasta) if line.startswith(">"))
        else:
            query_fasta = "/tmp/all_candidates_blast.fasta"
            n_seqs = combine_fastas(DATA_DIR, query_fasta)

    log.info(f"Query: {n_seqs} sequences, Mode: {args.mode}, E-value: {args.evalue}")

    output_tsv = str(OUTPUT_DIR / f"cas13_blast_{args.mode}.tsv")
    hits = run_blast(query_fasta, args.ref, output_tsv,
                     mode=args.mode, evalue=args.evalue)

    print("\n" + "=" * 100)
    print(f"{'COMPREHENSIVE CAS13 BLAST RESULTS':^100}")
    print(f"{'Reference DB: 8 proteins across Cas13a/b/d/bt subtypes':^100}")
    print("=" * 100)

    if not hits:
        print("\nNO HITS FOUND at e-value < {:.0e}".format(args.evalue))
        print("None of the candidates have detectable homology to any known Cas13 protein.")
        if args.mode == "blastp":
            print("\nTIP: Try blastx mode to search DNA contigs in all 6 reading frames:")
            print("  python blast_candidates.py --mode blastx --contigs /tmp/all_contigs.fasta")
    else:
        seen_queries = set()
        sorted_hits = sorted(hits, key=lambda x: x["evalue"])

        print(f"\n{'Query':45s} {'Subject':35s} {'%ID':>6s} {'Len':>5s} "
              f"{'E-value':>10s} {'Score':>7s} {'Coverage':>8s}")
        print("-" * 120)

        for h in sorted_hits:
            coverage = h["length"] / max(h["qlen"], 1) * 100
            seen_queries.add(h["query"])
            print(f"{h['query']:45s} {h['subject']:35s} {h['pident']:6.1f} "
                  f"{h['length']:5d} {h['evalue']:10.1e} {h['bitscore']:7.1f} "
                  f"{coverage:7.1f}%")

        print(f"\n{'=' * 100}")
        print(f"SUMMARY: {len(seen_queries)} sequences with Cas13 BLAST hits")
        print(f"Total hits: {len(sorted_hits)}")

        subtypes_hit = set(h["subject"].split("_")[0] for h in sorted_hits)
        print(f"Subtypes matched: {', '.join(sorted(subtypes_hit))}")

        strong = [q for q in seen_queries
                  if any(h["evalue"] < 1e-10 and h["pident"] > 20
                         for h in sorted_hits if h["query"] == q)]
        if strong:
            print(f"\nSTRONG Cas13 candidates (e<1e-10, >20% identity):")
            for q in sorted(strong):
                best = min([h for h in sorted_hits if h["query"] == q],
                           key=lambda x: x["evalue"])
                print(f"  {q} -> {best['subject']} ({best['pident']:.1f}% identity, "
                      f"e={best['evalue']:.1e})")

    print(f"\nFull results saved to: {output_tsv}")


if __name__ == "__main__":
    main()
