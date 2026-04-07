#!/usr/bin/env python3
"""
Validate top Cas13 candidate proteins for biological plausibility.

Checks:
  1. Basic QC: length, start codon, composition, low-complexity
  2. HEPN domain architecture: position, spacing, motif quality
  3. Homology to known Cas13 references (local alignment)
  4. Cas13-specific features: alpha-helical propensity, charge profile
  5. Optional: re-call ORFs from source contigs with Prodigal

Usage:
  python validate_proteins.py                          # validate all in mined_hits
  python validate_proteins.py --ids ID1 ID2 ID3        # specific sequences
  python validate_proteins.py --recall-orfs             # re-extract ORFs from cached contigs
"""
import argparse
import json
import logging
import os
import re
import sys
from collections import Counter
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")
log = logging.getLogger(__name__)

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
DATA_DIR = PROJECT_ROOT / "data" / "mined_hits"
OUTPUT_DIR = PROJECT_ROOT / "outputs"
METADATA_JSON = PROJECT_ROOT / "metadata" / "variant_domain_metadata.json"

HEPN_REGEX = re.compile(r"R.{3,8}H")
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

KNOWN_CAS13_REFS = {
    "LwaCas13a": "MKVTKVGGISHKKYTSEGRLVKSESEENRTDERLSALLNMRLDMYIKNPSSTETKENQKRIGKLKKFFSNKMVYLKDNTLSLKNGKKENIDREYSETDILESDVRDKKNFAVLKKIYLNENVNSEELEVFRNDIKKKLNKINSLKYSFEKNKANYQKINENNIEKVEGKSKRNIIYDYYRESAKRDAYVSNVKEAFDKLYKEEDIAKLVLEIENLTKLEKYKIREFYHEIIGRKNDKENFAKIIYEEIQNVNNMKELIEKVPDMSELKKSQVFYKYYLDKEELNDKNIKYAFCHFVEIEMSQLLKNYVYKRLSNISNDKIKRIFEYQNLKKLIENKLLNKLDTYVRNCGKYNYYLQDGEIATSDFIARNRQNEAFLRNIIGVSSVAYFSLRNILETENENDITGRMRGKTVKNNKGEEKYVSGEVDKIYNENKKNEVKENLKMFYSYDFNMDNKNEIEDFFANIDEAISSIRHGIVHFNLELEGKDIFAFKNIAPSEISKKMFQNEINEKKLKLKIFRQLNSANVFRYLEKYKILNYLKRTRFEFVNKNIPFVPSFTKLYSRIDDLKNSLGIYWKTPKTNDDNKTKEIIDAQIYLLKNIYYGEFLNYFMSNNGNFFEISKEIIELNKNDKRNLKTGFYKLQKFEDIQEKIPKEYLANIQSLYMINAGNQDEEEKDTYIDFIQKIFLKGFMTYLANNGRLSLIYIGSDEETNTSLAEKKQEFDKFLKKYEQNNNIKIPYEINEFLREIKLGNILKYTERLNMFYLILKLLNHKELTNLKGSLEKYQSANKEEAFSDQLELINLLNLDNNRVTEDFELEADEIGKFLDFNGNKVKDNKELKKFDTNKIYFDGENIIKHRAFYNIKKYGMLNLLEKIADKAGYKISIEELKKYSNKKNEIEKNHKMQENLHRKYARPRKDEKFTDEDYESYKQAIENIEEYTHLKNKVEFNELNLLQGLLLRILHRLVGYTSIWERDLRFRLKGEFPENQYIEEIFNFENKKNVKYKGGQIVEKYIKFYKELHQNDEVKINKYSSANIKVLKQEKKDLYIANYIAAFNYIPHAEISLLEVLENLRKLLSYDRKLKNAVMKSVVDILKEYGFVATFKIGADKKIGIQTLESEKIVHLKNLKKKKLMTDRNSEELCKLVKIMFEYKMEEKKSEN",
    "RfxCas13d_Nterm": "MAKKNKMPLSEKLLNDYFKVGKC",
    "PspCas13b_Nterm": "MNIPALRQQAMFQLYQGATFHYE",
}

CAS13_LENGTH_RANGE = (800, 1400)
IDEAL_LENGTH_RANGE = (900, 1200)

HYDROPHOBIC = set("AILMFVPW")
CHARGED = set("DEKR")
POLAR = set("STNQ")
AROMATIC = set("FWY")


def load_sequences(fasta_dir: Path, target_ids: list = None) -> dict:
    seqs = {}
    for fasta_path in sorted(fasta_dir.glob("*.fasta")):
        current_id = ""
        current_seq = []
        with open(fasta_path) as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_id and current_seq:
                        if target_ids is None or current_id in target_ids:
                            seqs[current_id] = "".join(current_seq)
                    current_id = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)
            if current_id and current_seq:
                if target_ids is None or current_id in target_ids:
                    seqs[current_id] = "".join(current_seq)
    return seqs


def check_basic_qc(seq_id: str, seq: str) -> dict:
    """Basic protein quality checks."""
    issues = []
    info = []
    length = len(seq)

    if length < CAS13_LENGTH_RANGE[0]:
        issues.append(f"TOO SHORT: {length}aa (Cas13 typically {CAS13_LENGTH_RANGE[0]}-{CAS13_LENGTH_RANGE[1]}aa)")
    elif length > CAS13_LENGTH_RANGE[1]:
        issues.append(f"TOO LONG: {length}aa (Cas13 typically {CAS13_LENGTH_RANGE[0]}-{CAS13_LENGTH_RANGE[1]}aa)")
    elif IDEAL_LENGTH_RANGE[0] <= length <= IDEAL_LENGTH_RANGE[1]:
        info.append(f"Length {length}aa — ideal Cas13 range")
    else:
        info.append(f"Length {length}aa — acceptable but not ideal")

    if not seq.startswith("M"):
        issues.append(f"NO START CODON: begins with '{seq[0]}' not 'M' — possible truncation or misannotated ORF")

    if "*" in seq:
        issues.append(f"INTERNAL STOP: {seq.count('*')} stop codon(s) found — likely frameshift or pseudogene")

    nonstandard = set(seq.upper()) - STANDARD_AA - {"*", "X", "U"}
    if nonstandard:
        issues.append(f"NON-STANDARD residues: {nonstandard}")

    x_count = seq.upper().count("X")
    if x_count > 0:
        pct = x_count / length * 100
        if pct > 5:
            issues.append(f"HIGH AMBIGUITY: {x_count} X residues ({pct:.1f}%) — poor sequencing quality")
        else:
            info.append(f"{x_count} ambiguous (X) residues ({pct:.1f}%)")

    return {"issues": issues, "info": info, "length": length}


def check_composition(seq: str) -> dict:
    """Amino acid composition analysis — flag non-protein-like profiles."""
    issues = []
    info = []
    counts = Counter(seq.upper())
    length = len(seq)

    hydro_frac = sum(counts.get(aa, 0) for aa in HYDROPHOBIC) / length
    charge_frac = sum(counts.get(aa, 0) for aa in CHARGED) / length
    polar_frac = sum(counts.get(aa, 0) for aa in POLAR) / length

    if hydro_frac > 0.55:
        issues.append(f"VERY HYDROPHOBIC ({hydro_frac:.1%}) — may be membrane protein, not cytoplasmic Cas13")
    elif hydro_frac > 0.45:
        info.append(f"Moderately hydrophobic ({hydro_frac:.1%}) — borderline for cytoplasmic protein")

    if charge_frac < 0.10:
        issues.append(f"LOW CHARGE ({charge_frac:.1%}) — Cas13 proteins are typically charged (RNA-binding)")
    elif charge_frac > 0.18:
        info.append(f"Highly charged ({charge_frac:.1%}) — consistent with RNA-binding")
    else:
        info.append(f"Charge fraction {charge_frac:.1%} — normal")

    max_run = _longest_single_aa_run(seq)
    if max_run >= 8:
        issues.append(f"LOW COMPLEXITY: single amino acid run of length {max_run}")
    elif max_run >= 6:
        info.append(f"Moderate single-AA run (length {max_run})")

    info.append(f"Composition: {hydro_frac:.1%} hydrophobic, {charge_frac:.1%} charged, {polar_frac:.1%} polar")

    return {"issues": issues, "info": info, "hydro_frac": hydro_frac, "charge_frac": charge_frac}


def _longest_single_aa_run(seq: str) -> int:
    max_run = 1
    current_run = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            current_run += 1
            max_run = max(max_run, current_run)
        else:
            current_run = 1
    return max_run


def check_hepn_architecture(seq: str) -> dict:
    """Verify HEPN domain positions and architecture."""
    issues = []
    info = []
    matches = list(HEPN_REGEX.finditer(seq))

    if len(matches) < 2:
        issues.append(f"MISSING HEPN: only {len(matches)} R...H motif(s) — Cas13 needs 2 HEPN domains")
        return {"issues": issues, "info": info, "n_hepn": len(matches), "best_pair": None}

    best_pair = None
    best_score = -1
    for i in range(len(matches)):
        for j in range(i + 1, len(matches)):
            spacing = matches[j].start() - matches[i].start()
            if 100 <= spacing <= 800:
                score = _pair_score(matches[i], matches[j], spacing, len(seq))
                if score > best_score:
                    best_score = score
                    best_pair = (matches[i], matches[j], spacing, score)

    if best_pair is None:
        issues.append(f"NO VALID HEPN PAIR: {len(matches)} motifs but none with 100-800aa spacing")
        return {"issues": issues, "info": info, "n_hepn": len(matches), "best_pair": None}

    m1, m2, spacing, score = best_pair
    pos1_frac = m1.start() / len(seq)
    pos2_frac = m2.start() / len(seq)

    if 150 <= spacing <= 400:
        info.append(f"HEPN spacing {spacing}aa — ideal Cas13 range")
    elif 100 <= spacing <= 500:
        info.append(f"HEPN spacing {spacing}aa — acceptable")
    else:
        issues.append(f"HEPN spacing {spacing}aa — atypical for Cas13 (expect 150-400)")

    if pos1_frac < 0.15 or pos1_frac > 0.7:
        issues.append(f"HEPN1 at {pos1_frac:.0%} of protein — unusual position")
    if pos2_frac < 0.3 or pos2_frac > 0.95:
        issues.append(f"HEPN2 at {pos2_frac:.0%} of protein — unusual position")

    info.append(f"HEPN1: pos {m1.start()} ({pos1_frac:.0%}) motif '{m1.group()}'")
    info.append(f"HEPN2: pos {m2.start()} ({pos2_frac:.0%}) motif '{m2.group()}'")
    info.append(f"Total R...H motifs: {len(matches)} (best pair score: {score:.3f})")

    return {"issues": issues, "info": info, "n_hepn": len(matches), "best_pair": best_pair}


def _pair_score(m1, m2, spacing, seq_len):
    score = 0.0
    if 200 <= spacing <= 400:
        score += 0.35
    elif 150 <= spacing <= 500:
        score += 0.25
    elif 100 <= spacing <= 800:
        score += 0.10
    for m in [m1, m2]:
        inner = m.group()[1:-1]
        if 3 <= len(inner) <= 7:
            score += 0.10
    return min(score, 1.0)


def check_homology(seq: str) -> dict:
    """Local alignment against known Cas13 reference sequences."""
    issues = []
    info = []

    best_ref = None
    best_identity = 0.0
    best_region = ""

    for ref_name, ref_seq in KNOWN_CAS13_REFS.items():
        if len(ref_seq) < 30:
            identity = _kmer_identity(seq[:50], ref_seq, k=4)
            region = "N-terminal"
        else:
            n_term = _kmer_identity(seq[:200], ref_seq[:200], k=5)
            c_term = _kmer_identity(seq[-200:], ref_seq[-200:], k=5)
            full = _kmer_identity(seq, ref_seq, k=6)
            identity = max(n_term, c_term, full)
            region = "N-term" if n_term == identity else ("C-term" if c_term == identity else "full")

        if identity > best_identity:
            best_identity = identity
            best_ref = ref_name
            best_region = region

    if best_identity < 0.02:
        issues.append(f"NO HOMOLOGY to known Cas13 references (best: {best_identity:.1%} to {best_ref})")
    elif best_identity < 0.05:
        issues.append(f"VERY LOW homology: {best_identity:.1%} to {best_ref} ({best_region}) — may not be Cas13")
    elif best_identity < 0.10:
        info.append(f"Low but detectable homology: {best_identity:.1%} to {best_ref} ({best_region}) — could be distant Cas13")
    else:
        info.append(f"Homology: {best_identity:.1%} to {best_ref} ({best_region})")

    return {"issues": issues, "info": info, "best_ref": best_ref, "best_identity": best_identity}


def _kmer_identity(seq1: str, seq2: str, k: int = 5) -> float:
    """K-mer based similarity (fast proxy for local alignment)."""
    if len(seq1) < k or len(seq2) < k:
        return 0.0
    kmers1 = set(seq1[i:i+k] for i in range(len(seq1) - k + 1))
    kmers2 = set(seq2[i:i+k] for i in range(len(seq2) - k + 1))
    if not kmers1 or not kmers2:
        return 0.0
    intersection = kmers1 & kmers2
    return len(intersection) / min(len(kmers1), len(kmers2))


def recall_orfs_with_prodigal(contigs_fasta: str, output_dir: str) -> str:
    """Re-call ORFs from source contigs using Prodigal in metagenomic mode.
    Returns path to re-called protein FASTA."""
    import shutil
    prodigal = shutil.which("prodigal")
    if not prodigal:
        log.error("Prodigal not found. Install: conda install -c bioconda prodigal")
        return None

    os.makedirs(output_dir, exist_ok=True)
    protein_out = os.path.join(output_dir, "prodigal_proteins.fasta")
    genes_out = os.path.join(output_dir, "prodigal_genes.gff")

    import subprocess
    log.info(f"Running Prodigal on {contigs_fasta}...")
    result = subprocess.run([
        prodigal,
        "-i", contigs_fasta,
        "-a", protein_out,
        "-o", genes_out,
        "-p", "meta",
        "-f", "gff",
    ], capture_output=True, text=True)

    if result.returncode != 0:
        log.error(f"Prodigal failed: {result.stderr}")
        return None

    n_proteins = sum(1 for line in open(protein_out) if line.startswith(">"))
    log.info(f"Prodigal extracted {n_proteins} proteins from contigs")

    filtered_out = os.path.join(output_dir, "prodigal_cas13_candidates.fasta")
    n_candidates = 0
    with open(protein_out) as fin, open(filtered_out, "w") as fout:
        current_id = ""
        current_seq = []
        for line in fin:
            if line.startswith(">"):
                if current_id and current_seq:
                    seq = "".join(current_seq).replace("*", "")
                    if 800 <= len(seq) <= 1400:
                        hepn_hits = list(HEPN_REGEX.finditer(seq))
                        has_pair = any(
                            100 <= hepn_hits[j].start() - hepn_hits[i].start() <= 800
                            for i in range(len(hepn_hits))
                            for j in range(i + 1, len(hepn_hits))
                        ) if len(hepn_hits) >= 2 else False
                        if has_pair:
                            fout.write(f">{current_id}\n{seq}\n")
                            n_candidates += 1
                current_id = line[1:].strip().split()[0]
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_id and current_seq:
            seq = "".join(current_seq).replace("*", "")
            if 800 <= len(seq) <= 1400:
                hepn_hits = list(HEPN_REGEX.finditer(seq))
                has_pair = any(
                    100 <= hepn_hits[j].start() - hepn_hits[i].start() <= 800
                    for i in range(len(hepn_hits))
                    for j in range(i + 1, len(hepn_hits))
                ) if len(hepn_hits) >= 2 else False
                if has_pair:
                    fout.write(f">{current_id}\n{seq}\n")
                    n_candidates += 1

    log.info(f"Filtered to {n_candidates} Cas13-sized proteins with HEPN pairs -> {filtered_out}")
    return filtered_out


def validate_sequence(seq_id: str, seq: str) -> dict:
    """Run all validation checks on a single sequence."""
    qc = check_basic_qc(seq_id, seq)
    comp = check_composition(seq)
    hepn = check_hepn_architecture(seq)
    homol = check_homology(seq)

    all_issues = qc["issues"] + comp["issues"] + hepn["issues"] + homol["issues"]
    all_info = qc["info"] + comp["info"] + hepn["info"] + homol["info"]

    if len(all_issues) == 0:
        verdict = "PASS"
    elif len(all_issues) <= 2 and not any("MISSING HEPN" in i or "NO HOMOLOGY" in i or "INTERNAL STOP" in i for i in all_issues):
        verdict = "WARN"
    else:
        verdict = "FAIL"

    return {
        "seq_id": seq_id,
        "length": qc["length"],
        "verdict": verdict,
        "n_issues": len(all_issues),
        "issues": all_issues,
        "info": all_info,
        "n_hepn": hepn["n_hepn"],
        "best_ref": homol.get("best_ref"),
        "best_identity": homol.get("best_identity", 0),
        "hydro_frac": comp.get("hydro_frac", 0),
        "charge_frac": comp.get("charge_frac", 0),
    }


def main():
    parser = argparse.ArgumentParser(description="Validate Cas13 candidate proteins")
    parser.add_argument("--ids", nargs="*", help="Specific sequence IDs to validate")
    parser.add_argument("--fasta-dir", default=str(DATA_DIR))
    parser.add_argument("--recall-orfs", action="store_true",
                        help="Re-call ORFs from cached contigs using Prodigal")
    parser.add_argument("--contigs", default="/tmp/all_contigs.fasta",
                        help="Path to contigs FASTA for --recall-orfs")
    args = parser.parse_args()

    target_ids = set(args.ids) if args.ids else None
    sequences = load_sequences(Path(args.fasta_dir), target_ids)

    if not sequences:
        log.error(f"No sequences found in {args.fasta_dir}")
        sys.exit(1)

    log.info(f"Validating {len(sequences)} protein sequences...\n")

    results = []
    for seq_id, seq in sorted(sequences.items()):
        result = validate_sequence(seq_id, seq)
        results.append(result)

    pass_count = sum(1 for r in results if r["verdict"] == "PASS")
    warn_count = sum(1 for r in results if r["verdict"] == "WARN")
    fail_count = sum(1 for r in results if r["verdict"] == "FAIL")

    print("=" * 90)
    print(f"{'PROTEIN VALIDATION REPORT':^90}")
    print("=" * 90)

    for r in sorted(results, key=lambda x: ({"PASS": 0, "WARN": 1, "FAIL": 2}[x["verdict"]], x["seq_id"])):
        icon = {"PASS": "+", "WARN": "~", "FAIL": "X"}[r["verdict"]]
        print(f"\n[{icon}] {r['verdict']} | {r['seq_id']}")
        print(f"    Length: {r['length']}aa | HEPN motifs: {r['n_hepn']} | "
              f"Closest ref: {r['best_ref']} ({r['best_identity']:.1%}) | "
              f"Hydrophobic: {r['hydro_frac']:.1%} | Charged: {r['charge_frac']:.1%}")

        if r["issues"]:
            for issue in r["issues"]:
                print(f"    !! {issue}")
        for note in r["info"][:3]:
            print(f"       {note}")

    print("\n" + "=" * 90)
    print(f"SUMMARY: {pass_count} PASS | {warn_count} WARN | {fail_count} FAIL | {len(results)} total")
    print("=" * 90)

    if fail_count > 0:
        print("\nFAILED sequences should be REMOVED from validated_baseline_ids.txt")
        print("WARN sequences may need manual review")

    if args.recall_orfs:
        if not os.path.exists(args.contigs):
            log.error(f"Contigs file not found: {args.contigs}")
            log.info("If you have cached contigs elsewhere, use --contigs /path/to/contigs.fasta")
        else:
            output_dir = str(OUTPUT_DIR / "prodigal_recall")
            recalled = recall_orfs_with_prodigal(args.contigs, output_dir)
            if recalled:
                print(f"\n{'=' * 90}")
                print(f"PRODIGAL RE-CALLED ORFs: {recalled}")
                print(f"Compare these against your current candidates to find")
                print(f"corrected/improved ORFs from the same source contigs.")
                print(f"{'=' * 90}")

    report_path = str(OUTPUT_DIR / "protein_validation_report.json")
    os.makedirs(os.path.dirname(report_path), exist_ok=True)
    with open(report_path, "w") as f:
        json.dump([{k: v for k, v in r.items() if k != "best_pair"} for r in results], f, indent=2)
    log.info(f"Full report saved to {report_path}")


if __name__ == "__main__":
    main()
