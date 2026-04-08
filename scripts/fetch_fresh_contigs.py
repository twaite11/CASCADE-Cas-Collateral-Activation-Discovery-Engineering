#!/usr/bin/env python3
"""
Fetch the newest bacterial and archaeal WGS contigs from NCBI for Cas13 mining.

Strategy: Query NCBI Assembly database for the most recently submitted
bacterial/archaeal genome assemblies, download their contig-level nucleotide
sequences, and output a combined FASTA ready for mining_v2.py.

Usage:
  python fetch_fresh_contigs.py --count 4000 --output /workspace/CASCADE/data/fresh_contigs.fasta
  python fetch_fresh_contigs.py --count 4000 --output fresh.fasta --run-mining

Requires: biopython (pip install biopython)
"""
import argparse
import logging
import os
import sys
import time
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent


def fetch_fresh_assemblies(count: int = 4000, email: str = "founder@senarybio.com",
                           min_contig_len: int = 50000) -> list[str]:
    """
    Query NCBI Assembly DB for the newest bacterial + archaeal WGS assemblies.
    Returns list of WGS nucleotide accessions (NZ_* contig accessions).
    """
    from Bio import Entrez
    Entrez.email = email
    Entrez.api_key = os.environ.get("NCBI_API_KEY", None)

    rate_delay = 0.11 if Entrez.api_key else 0.35

    target_per_kingdom = count // 2
    all_accessions = []

    for kingdom, term_fragment in [
        ("Bacteria", '"Bacteria"[Organism]'),
        ("Archaea", '"Archaea"[Organism]'),
    ]:
        log.info(f"Searching NCBI Assembly for newest {kingdom} genomes...")

        search_term = (
            f'{term_fragment} AND "latest"[filter] AND '
            f'"contig level"[filter] OR "scaffold level"[filter] OR "complete genome"[filter]'
        )

        handle = Entrez.esearch(
            db="assembly",
            term=search_term,
            retmax=target_per_kingdom,
            sort="Date Released",  # newest first
            usehistory="y",
        )
        search_results = Entrez.read(handle)
        handle.close()
        time.sleep(rate_delay)

        id_list = search_results.get("IdList", [])
        total_found = int(search_results.get("Count", 0))
        log.info(f"  {kingdom}: {total_found} total assemblies, fetching {len(id_list)} newest")

        if not id_list:
            continue

        batch_size = 200
        for start in range(0, len(id_list), batch_size):
            batch = id_list[start:start + batch_size]
            log.info(f"  Fetching assembly metadata batch {start // batch_size + 1} "
                     f"({len(batch)} assemblies)...")

            try:
                handle = Entrez.esummary(db="assembly", id=",".join(batch), retmax=batch_size)
                summaries = Entrez.read(handle)
                handle.close()
                time.sleep(rate_delay)
            except Exception as e:
                log.warning(f"  Failed to fetch batch: {e}")
                continue

            doc_sums = summaries.get("DocumentSummarySet", {}).get("DocumentSummary", [])
            for doc in doc_sums:
                accession = doc.get("AssemblyAccession", "")
                wgs_project = doc.get("WGS", "")
                gbrs_paired = doc.get("GbrsAccession", "")
                refseq_acc = doc.get("RsUid", "")
                synonym = doc.get("Synonym", {})
                genbank_acc = synonym.get("Genbank", "") if isinstance(synonym, dict) else ""
                refseq_ftp = doc.get("FtpPath_RefSeq", "")
                genbank_ftp = doc.get("FtpPath_GenBank", "")

                chosen_acc = genbank_acc or accession
                if chosen_acc:
                    all_accessions.append(chosen_acc)

            log.info(f"  Running total: {len(all_accessions)} assembly accessions")

    log.info(f"Collected {len(all_accessions)} assembly accessions total")
    return all_accessions[:count]


def download_contigs_via_efetch(accessions: list[str], output_fasta: str,
                                 email: str = "founder@senarybio.com",
                                 batch_size: int = 50, min_len: int = 20000) -> int:
    """
    Download nucleotide contigs for each assembly accession.
    Uses Entrez elink (assembly -> nucleotide) then efetch in batches.
    Returns total contig count written.
    """
    from Bio import Entrez, SeqIO
    Entrez.email = email
    Entrez.api_key = os.environ.get("NCBI_API_KEY", None)

    rate_delay = 0.11 if Entrez.api_key else 0.35
    total_written = 0
    total_bases = 0

    os.makedirs(os.path.dirname(os.path.abspath(output_fasta)), exist_ok=True)

    with open(output_fasta, "w") as out_f:
        for idx, acc in enumerate(accessions):
            if total_written >= 4500:
                log.info(f"Reached target contig count ({total_written}), stopping fetch")
                break

            if idx % 100 == 0 and idx > 0:
                log.info(f"Progress: {idx}/{len(accessions)} assemblies processed, "
                         f"{total_written} contigs written, {total_bases / 1e9:.2f} Gbp")

            try:
                handle = Entrez.esearch(
                    db="nucleotide",
                    term=f"{acc}[Assembly]",
                    retmax=50,
                )
                nuc_results = Entrez.read(handle)
                handle.close()
                time.sleep(rate_delay)

                nuc_ids = nuc_results.get("IdList", [])
                if not nuc_ids:
                    continue

                handle = Entrez.efetch(
                    db="nucleotide",
                    id=",".join(nuc_ids[:20]),
                    rettype="fasta",
                    retmode="text",
                )
                for rec in SeqIO.parse(handle, "fasta"):
                    seq_len = len(rec.seq)
                    if seq_len >= min_len:
                        out_f.write(f">{rec.id}\n{str(rec.seq)}\n")
                        total_written += 1
                        total_bases += seq_len
                handle.close()
                time.sleep(rate_delay)

            except Exception as e:
                log.warning(f"  Failed on {acc}: {e}")
                time.sleep(1)
                continue

    log.info(f"Download complete: {total_written} contigs, {total_bases / 1e9:.2f} Gbp")
    return total_written


def download_contigs_via_datasets(accessions: list[str], output_fasta: str,
                                   batch_size: int = 100) -> int:
    """
    Use NCBI datasets CLI (much faster) if available.
    Falls back to efetch if datasets not installed.
    """
    import shutil
    import subprocess
    import tempfile
    import zipfile

    datasets_bin = shutil.which("datasets")
    if not datasets_bin:
        log.info("NCBI datasets CLI not found, falling back to Entrez efetch (slower)")
        return -1

    log.info(f"Using NCBI datasets CLI for fast bulk download")
    os.makedirs(os.path.dirname(os.path.abspath(output_fasta)), exist_ok=True)

    total_written = 0
    total_bases = 0

    with open(output_fasta, "w") as out_f:
        for start in range(0, len(accessions), batch_size):
            batch = accessions[start:start + batch_size]
            log.info(f"Batch {start // batch_size + 1}: downloading {len(batch)} assemblies...")

            with tempfile.TemporaryDirectory() as tmpdir:
                acc_file = os.path.join(tmpdir, "accs.txt")
                with open(acc_file, "w") as f:
                    f.write("\n".join(batch))

                zip_path = os.path.join(tmpdir, "dataset.zip")
                try:
                    result = subprocess.run(
                        [datasets_bin, "download", "genome", "accession",
                         "--inputfile", acc_file,
                         "--include", "genome",
                         "--filename", zip_path],
                        capture_output=True, text=True, timeout=600,
                    )
                    if result.returncode != 0:
                        log.warning(f"datasets CLI error: {result.stderr[:200]}")
                        continue
                except subprocess.TimeoutExpired:
                    log.warning(f"Batch timed out, skipping")
                    continue

                if not os.path.exists(zip_path):
                    continue

                try:
                    with zipfile.ZipFile(zip_path) as zf:
                        for name in zf.namelist():
                            if name.endswith(".fna") or name.endswith(".fasta"):
                                from Bio import SeqIO
                                import io
                                with zf.open(name) as fna:
                                    text = io.TextIOWrapper(fna)
                                    for rec in SeqIO.parse(text, "fasta"):
                                        if len(rec.seq) >= 20000:
                                            out_f.write(f">{rec.id}\n{str(rec.seq)}\n")
                                            total_written += 1
                                            total_bases += len(rec.seq)
                except Exception as e:
                    log.warning(f"Failed to parse zip: {e}")

            log.info(f"  Running total: {total_written} contigs, {total_bases / 1e9:.2f} Gbp")

            if total_written >= 4500:
                log.info("Reached target contig count, stopping")
                break

    log.info(f"Download complete: {total_written} contigs, {total_bases / 1e9:.2f} Gbp")
    return total_written


def main():
    parser = argparse.ArgumentParser(
        description="Fetch newest bacterial/archaeal contigs from NCBI for Cas13 mining"
    )
    parser.add_argument("--count", type=int, default=4000,
                        help="Target number of assemblies to query (default: 4000)")
    parser.add_argument("--output", type=str,
                        default=str(PROJECT_ROOT / "data" / "fresh_contigs.fasta"),
                        help="Output FASTA path")
    parser.add_argument("--email", type=str, default="founder@senarybio.com")
    parser.add_argument("--min-contig-len", type=int, default=20000,
                        help="Minimum contig length to keep (bp, default: 20000)")
    parser.add_argument("--run-mining", action="store_true",
                        help="Automatically run mining_v2.py after download")
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    log.info("=" * 70)
    log.info("CASCADE Fresh Contig Fetcher — Mining the Newest Genomes")
    log.info("=" * 70)
    log.info(f"Target assemblies: {args.count}")
    log.info(f"Output: {args.output}")

    assembly_accs = fetch_fresh_assemblies(
        count=args.count,
        email=args.email,
    )

    if not assembly_accs:
        log.error("No assemblies found. Check network / NCBI availability.")
        sys.exit(1)

    n_written = download_contigs_via_datasets(assembly_accs, args.output)
    if n_written <= 0:
        log.info("datasets CLI produced no contigs, falling back to Entrez efetch...")
        n_written = download_contigs_via_efetch(
            assembly_accs, args.output,
            email=args.email,
            min_len=args.min_contig_len,
        )

    if n_written == 0:
        log.error("No contigs downloaded. Check accessions / network.")
        sys.exit(1)

    log.info(f"\nFinal: {n_written} contigs saved to {args.output}")

    if args.run_mining:
        log.info("\n" + "=" * 70)
        log.info("Starting mining_v2.py on fresh contigs...")
        log.info("=" * 70)
        mining_script = SCRIPT_DIR / "mining_v2.py"
        mining_output = str(PROJECT_ROOT / "outputs" / "mining_v2_fresh")
        os.makedirs(mining_output, exist_ok=True)

        import subprocess
        cmd = [
            sys.executable, str(mining_script),
            "--contigs", args.output,
            "--output-dir", mining_output,
            "--threads", str(args.threads),
        ]
        log.info(f"Running: {' '.join(cmd)}")
        result = subprocess.run(cmd)
        sys.exit(result.returncode)
    else:
        log.info(f"\nTo run mining:")
        log.info(f"  python mining_v2.py --contigs {args.output} "
                 f"--output-dir ../outputs/mining_v2_fresh --threads {args.threads}")


if __name__ == "__main__":
    main()
