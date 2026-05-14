"""Fetch a list of NCBI nucleotide accessions to a single FASTA, with cache.

Used for the 2026-05-13 re-mine of Campaign 1 (Listeria booriae) and
Campaign 2 (Bacteroides + Flavobacterium + Leptotrichia) contigs.
"""
from __future__ import annotations

import argparse
import logging
import time
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "data" / "contig_cache"


def fetch_one(acc: str, email: str) -> str | None:
    CACHE.mkdir(parents=True, exist_ok=True)
    cache_file = CACHE / f"{acc}.fasta"
    if cache_file.exists() and cache_file.stat().st_size > 1000:
        log.info(f"  {acc}: cache hit ({cache_file.stat().st_size:,} bytes)")
        return cache_file.read_text()
    try:
        from Bio import Entrez, SeqIO
        Entrez.email = email
        handle = Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text")
        records = list(SeqIO.parse(handle, "fasta"))
        handle.close()
        if not records:
            log.warning(f"  {acc}: no records returned")
            return None
        text = "".join(f">{r.id}\n{str(r.seq)}\n" for r in records)
        cache_file.write_text(text)
        bp = sum(len(r.seq) for r in records)
        log.info(f"  {acc}: fetched {len(records)} record(s), {bp:,} bp")
        time.sleep(0.5)
        return text
    except Exception as e:
        log.error(f"  {acc}: fetch failed: {e}")
        return None


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--accessions", nargs="+", required=True,
                   help="NCBI nucleotide accessions to fetch")
    p.add_argument("--out", required=True, type=Path,
                   help="Output combined FASTA path")
    p.add_argument("--email", default="founder@senarybio.com")
    args = p.parse_args()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    n_ok = 0
    with args.out.open("w") as fh:
        for acc in args.accessions:
            text = fetch_one(acc, args.email)
            if text:
                fh.write(text)
                n_ok += 1
    bp = sum(len(line) for line in args.out.read_text().splitlines() if not line.startswith(">"))
    log.info(f"Wrote {n_ok}/{len(args.accessions)} contigs to {args.out} ({bp:,} bp)")


if __name__ == "__main__":
    main()
