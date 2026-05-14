"""Quarantine all non-Cas13 data from the active CASCADE pipeline.

Moves files into quarantine/non_cas13_audit_<STAMP>/ with full audit trail.
Trims FASTAs/CSVs in place (with .bak originals saved alongside).
Wipes the SQLite catalog (originals backed up).

Run from repo root with the project's .venv python:
    .venv\\Scripts\\python.exe scripts\\quarantine_non_cas13.py
"""
from __future__ import annotations

import csv
import shutil
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STAMP = "20260513"
QDIR = ROOT / "quarantine" / f"non_cas13_audit_{STAMP}"
QDIR.mkdir(parents=True, exist_ok=True)

LOG: list[str] = []


def log(msg: str) -> None:
    print(msg)
    LOG.append(msg)


def quarantine_file(src: Path, reason: str) -> None:
    if not src.exists():
        log(f"  SKIP (not present): {src}")
        return
    rel = src.relative_to(ROOT)
    dst = QDIR / rel
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.is_dir():
        shutil.move(str(src), str(dst))
        log(f"  MOVE dir : {rel}  ->  {dst.relative_to(ROOT)}   [{reason}]")
    else:
        shutil.move(str(src), str(dst))
        log(f"  MOVE file: {rel}  ->  {dst.relative_to(ROOT)}   [{reason}]")


# Verified-real Cas13a IDs (KEEP)
CONFIRMED_CAS13A = {
    "NZ_JAASWF010000012.1_ORF_f1_6252_HIGH",
    "NZ_JAAROR010000001.1_ORF_f1_65475_HIGH",
    "NZ_JAARYF010000017.1_ORF_f1_4771_HIGH",
}


def trim_fasta(path: Path, keep_ids: set[str]) -> None:
    if not path.exists():
        log(f"  SKIP (not present): {path}")
        return
    bak = path.with_suffix(path.suffix + ".bak_audit_" + STAMP)
    shutil.copy2(path, bak)
    log(f"  BACKUP: {path.relative_to(ROOT)}  ->  {bak.relative_to(ROOT)}")

    kept: list[tuple[str, list[str]]] = []
    cur: tuple[str, list[str]] | None = None
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    kept.append(cur)
                cur = (line.rstrip()[1:].split()[0], [])
            elif cur is not None and line.strip():
                cur[1].append(line.rstrip())
        if cur is not None:
            kept.append(cur)

    surviving = [(h, s) for (h, s) in kept if h in keep_ids]
    dropped = [h for (h, _) in kept if h not in keep_ids]
    with path.open("w") as fh:
        for h, s in surviving:
            fh.write(f">{h}\n")
            fh.write("\n".join(s) + "\n")
    log(f"  TRIM fasta {path.relative_to(ROOT)}  kept={len(surviving)}  dropped={len(dropped)}")
    for d in dropped:
        log(f"    dropped: {d}")


def trim_csv(path: Path, keep_ids: set[str], id_col: str) -> None:
    if not path.exists():
        log(f"  SKIP (not present): {path}")
        return
    bak = path.with_suffix(path.suffix + ".bak_audit_" + STAMP)
    shutil.copy2(path, bak)
    log(f"  BACKUP: {path.relative_to(ROOT)}  ->  {bak.relative_to(ROOT)}")

    with path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)
        fns = reader.fieldnames or []

    surviving = [r for r in rows if r.get(id_col) in keep_ids]
    dropped = [r.get(id_col) for r in rows if r.get(id_col) not in keep_ids]
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fns)
        w.writeheader()
        w.writerows(surviving)
    log(f"  TRIM csv {path.relative_to(ROOT)}  kept={len(surviving)}  dropped={len(dropped)}")
    for d in dropped:
        log(f"    dropped: {d}")


def main() -> None:
    log(f"# CASCADE non-Cas13 data quarantine -- {datetime.now(timezone.utc).isoformat()}")
    log(f"# Quarantine dir: {QDIR.relative_to(ROOT)}")
    log("")

    # 1. Trim curated mined_hits set
    log("## 1. Trim data/mined_hits/confirmed_novel_cas13.* to confirmed Cas13a only")
    trim_fasta(ROOT / "data" / "mined_hits" / "confirmed_novel_cas13.fasta", CONFIRMED_CAS13A)
    trim_csv(ROOT / "data" / "mined_hits" / "confirmed_novel_cas13_metadata.csv",
             CONFIRMED_CAS13A, "sequence_id")
    log("")

    # 2. Move Campaign 2 fresh-mining raw outputs to data/recoverable_inputs/campaign2
    log("## 2. Move Campaign 2 fresh-mining raw outputs to data/recoverable_inputs/campaign2/")
    log("##    (raw mining_v2_hits.fasta, diamond_hits.tsv, minced_output.gff PRESERVED)")
    camp2 = ROOT / "outputs" / "mining_v2_fresh_v2"
    if camp2.exists():
        keep_dir = ROOT / "data" / "recoverable_inputs" / "campaign2"
        keep_dir.mkdir(parents=True, exist_ok=True)
        for f in list(camp2.iterdir()):
            if f.is_file():
                shutil.move(str(f), str(keep_dir / f.name))
                log(f"  PRESERVE: {f.relative_to(ROOT)} -> {(keep_dir / f.name).relative_to(ROOT)}")
        try:
            camp2.rmdir()
            log(f"  Removed empty {camp2.relative_to(ROOT)}")
        except OSError:
            pass
    log("")

    # 3. Quarantine Campaign 1 numeric-ID mining outputs
    log("## 3. Quarantine Campaign 1 numeric-ID mining_v2/ outputs")
    camp1 = ROOT / "outputs" / "mining_v2"
    if camp1.exists():
        quarantine_file(camp1, "Campaign 1 numeric-ID outputs (110 catalog entries, ~58% non-Cas13)")
    log("")

    # 4. JSONs: keep only those for confirmed Cas13a + positive control
    log("## 4. Quarantine eval JSONs not for confirmed Cas13a")
    jsons = ROOT / "jsons"
    if jsons.exists():
        keep_prefixes = tuple(s.replace("_HIGH", "") for s in CONFIRMED_CAS13A)
        for f in sorted(jsons.iterdir()):
            if not f.is_file():
                continue
            keep = (
                any(f.name.startswith(p) for p in keep_prefixes)
                or f.name == "Cas13a_positive_control.json"
            )
            if not keep:
                quarantine_file(f, "non-confirmed-Cas13a eval JSON")
    log("")

    # 5. Phase1 screening structures: keep only confirmed Cas13a
    log("## 5. Quarantine phase1_screening structures not for confirmed Cas13a")
    p1 = ROOT / "outputs" / "phase1_screening"
    if p1.exists():
        keep_prefixes = tuple(s.replace("_HIGH", "") for s in CONFIRMED_CAS13A)
        for f in sorted(p1.iterdir()):
            if f.is_file():
                keep = any(f.name.startswith(p) for p in keep_prefixes)
                if not keep:
                    quarantine_file(f, "non-confirmed-Cas13a phase1 structure")
    log("")

    # 6. Optimized switches with unknown lineage origin
    log("## 6. Quarantine optimized_switches whose baseline cannot be tied to confirmed Cas13a")
    opt = ROOT / "outputs" / "optimized_switches"
    if opt.exists():
        for f in sorted(opt.iterdir()):
            if f.is_file():
                name = f.name
                keep = (
                    name.startswith("NZ_JAASWF")
                    or name.startswith("NZ_JAAROR")
                    or name.startswith("NZ_JAARYF")
                )
                if not keep:
                    quarantine_file(f, "evolved variant from unverified lineage (L*/3174*)")
            elif f.is_dir() and f.name == ".ipynb_checkpoints":
                quarantine_file(f, "ipynb checkpoints noise")
    log("")

    # 7. Quarantine deep_hits_*_metadata.csv (Campaign 1 deep-mining metadata)
    log("## 7. Quarantine deep_hits_* Campaign 1 deep-mining metadata CSVs")
    mh = ROOT / "data" / "mined_hits"
    deep_csvs = sorted(mh.glob("deep_hits_*_metadata.csv"))
    log(f"  found {len(deep_csvs)} deep_hits_*_metadata.csv files")
    for f in deep_csvs:
        quarantine_file(f, "Campaign 1 deep-mining (numeric IDs, mostly non-Cas13)")
    deep_fastas = sorted(mh.glob("deep_hits_*.fasta"))
    for f in deep_fastas:
        quarantine_file(f, "Campaign 1 deep-mining FASTA (numeric IDs)")
    log("")

    # 8. SQLite catalog: back up + wipe contents (preserve schema)
    log("## 8. Wipe metadata/cas13_variants.db (preserve schema, archive original)")
    db = ROOT / "metadata" / "cas13_variants.db"
    shm = ROOT / "metadata" / "cas13_variants.db-shm"
    wal = ROOT / "metadata" / "cas13_variants.db-wal"
    if db.exists():
        bak_dir = QDIR / "metadata"
        bak_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(db, bak_dir / f"cas13_variants.db.bak_audit_{STAMP}")
        if shm.exists():
            shutil.copy2(shm, bak_dir / shm.name)
        if wal.exists():
            shutil.copy2(wal, bak_dir / wal.name)
        log(f"  BACKUP DB to {bak_dir.relative_to(ROOT)}/cas13_variants.db.bak_audit_{STAMP}")
        con = sqlite3.connect(db)
        n_before = con.execute("SELECT COUNT(*) FROM variants").fetchone()[0]
        con.execute("DELETE FROM variants")
        con.commit()
        con.execute("VACUUM")
        n_after = con.execute("SELECT COUNT(*) FROM variants").fetchone()[0]
        con.close()
        for f in (shm, wal):
            if f.exists():
                f.unlink()
                log(f"  removed stale {f.relative_to(ROOT)}")
        log(f"  WIPED variants table: {n_before} -> {n_after} rows")
    log("")

    # 9. Quarantine variant_domain_metadata.json (it indexes the wiped variants)
    log("## 9. Quarantine metadata/variant_domain_metadata.json (indexes wiped catalog)")
    vdm = ROOT / "metadata" / "variant_domain_metadata.json"
    if vdm.exists():
        quarantine_file(vdm, "indexes wiped Campaign 1 catalog; will be regenerated")
    log("")

    # 10. Drop the audit verification scripts created in this session
    log("## 10. Move temporary audit scripts to quarantine for reference")
    for scratch in ("verify_candidates.py", "inspect_db.py", "audit_variants_db.py",
                    "audit_variants_classification.csv"):
        p = ROOT / scratch
        if p.exists():
            quarantine_file(p, "temporary audit script (kept for traceability)")
    log("")

    # Write the manifest
    manifest = QDIR / "MANIFEST.md"
    confirmed = "\n".join(f"- {x}" for x in sorted(CONFIRMED_CAS13A))
    log_text = "\n".join(LOG)
    manifest_text = (
        f"# CASCADE non-Cas13 quarantine manifest -- {STAMP}\n\n"
        f"Generated: {datetime.now(timezone.utc).isoformat()}\n\n"
        f"## Why\n\n"
        f"Programmatic re-audit of `CASCADE_CANDIDATE_REPORT.md` (April 8, 2026)\n"
        f"and the entire `metadata/cas13_variants.db` catalog (110 entries) showed\n"
        f"that the majority of pipeline-classified Cas13 hits are mis-annotated\n"
        f"non-Cas13 enzymes. Diagnostic motifs prove these are: glycoside\n"
        f"hydrolase family 3 (KHFPGHGD), cobalamin-dependent methionine synthase\n"
        f"(DGAMGTM + MTGMDEVGRLF), metallo-beta-lactamase / RNase Z (LTHLHSDHV),\n"
        f"AAA+ ATPases (Walker A), sugar isomerases, secreted lipoproteins\n"
        f"(lipobox), and one Cas9 (LGLDLGTNSIGW). The pipeline bugs that caused\n"
        f"these false positives are tracked as audit findings B-2, B-10, B-13,\n"
        f"B-14 (loose `R.{{3,8}}H` HEPN regex), and B-9.\n\n"
        f"## Confirmed real Cas13a (PRESERVED in active pipeline)\n\n"
        f"{confirmed}\n\n"
        f"All three are *Listeria booriae* Cas13a (1052-1064 aa), 78-81%% to\n"
        f"patented Cas13a (US 10,266,886 Seq 11), 5+ canonical R-X(4-6)-H HEPN\n"
        f"motifs, no diagnostic non-Cas13 fold signatures.\n\n"
        f"## Action log\n\n"
        f"```\n{log_text}\n```\n"
    )
    manifest.write_text(manifest_text, encoding="utf-8")
    log(f"\nManifest written: {manifest.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
