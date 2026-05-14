"""Build ``metadata/variant_domain_metadata.json`` for the three confirmed
novel Cas13 enzymes so the active-learning orchestrator can pick them up.

The orchestrator expects a flat dict keyed by ``sequence_id``:

  {
    "<sequence_id>": {
      "sequence_length": <int>,
      "domains": {
        "HEPN1": {"start": <int>, "end": <int>},
        "HEPN2": {"start": <int>, "end": <int>}
      },
      "crRNA_repeat_used": "<RNA sequence, U-form>",
      "subtype": "cas13a"
    },
    ...
  }

Inputs (committed in-repo):
  data/mined_hits/confirmed_novel_cas13.fasta       -- protein sequences
  data/mined_hits/confirmed_novel_cas13_metadata.csv -- DR sequences

Window convention matches build_metadata_override_for_evolved in
scripts/evolution_orchestrator.py (HEPN1 = [center-30, center+80],
HEPN2 = [max(prev_end, center-30), center+80]).
"""
from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from evolution_orchestrator import _select_hepn_pair  # noqa: E402

FASTA = ROOT / "data" / "mined_hits" / "confirmed_novel_cas13.fasta"
CSV = ROOT / "data" / "mined_hits" / "confirmed_novel_cas13_metadata.csv"
OUT = ROOT / "metadata" / "variant_domain_metadata.json"

HEPN_MOTIF = re.compile(r"R.{4,6}H")  # strict Cas13 HEPN (B-3)


def parse_fasta(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    cur_name: str | None = None
    cur_lines: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        if raw.startswith(">"):
            if cur_name is not None:
                out[cur_name] = "".join(cur_lines)
            cur_name = raw[1:].split()[0]
            cur_lines = []
        else:
            cur_lines.append(raw.strip())
    if cur_name is not None:
        out[cur_name] = "".join(cur_lines)
    return out


def build() -> dict[str, dict]:
    seqs = parse_fasta(FASTA)
    dr_map: dict[str, str] = {}
    with CSV.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            dr_map[row["sequence_id"]] = row["repeat_domains"].strip()

    out: dict[str, dict] = {}
    for sid, prot in seqs.items():
        if sid not in dr_map:
            print(f"WARN  {sid}: no DR in metadata csv -- skipping")
            continue
        pair = _select_hepn_pair(prot, HEPN_MOTIF)
        if pair is None:
            print(f"WARN  {sid}: no valid HEPN pair found -- skipping")
            continue
        h1c = pair[0].start()
        h2c = pair[1].start()
        h1_start = max(0, h1c - 30)
        h1_end = h1c + 80
        h2_start = max(h1_end, h2c - 30)
        h2_end = h2c + 80
        dr_dna = dr_map[sid].upper()
        dr_rna = dr_dna.replace("T", "U")
        out[sid] = {
            "sequence_length": len(prot),
            "domains": {
                "HEPN1": {"start": h1_start, "end": h1_end},
                "HEPN2": {"start": h2_start, "end": h2_end},
            },
            "crRNA_repeat_used": dr_rna,
            "subtype": "cas13a",
        }
        print(
            f"OK    {sid:50s} len={len(prot):4d}  "
            f"HEPN1=R{pair[0].start()+1}-H{pair[0].end()}  "
            f"HEPN2=R{pair[1].start()+1}-H{pair[1].end()}  "
            f"sep={pair[1].start() - pair[0].start()}  "
            f"DR_nt={len(dr_rna)}"
        )
    return out


def main() -> int:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    data = build()
    if not data:
        print("ERROR no enzymes produced -- not writing metadata file", file=sys.stderr)
        return 1
    OUT.write_text(json.dumps(data, indent=2), encoding="utf-8")
    print(f"\nWrote {OUT} with {len(data)} enzymes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
