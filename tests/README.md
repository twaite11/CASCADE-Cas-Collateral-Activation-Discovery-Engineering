# 🧪 CASCADE Tests

Unit tests for all pipeline components. **No GPU, Protenix, or PXDesign installation required.**

---

## Setup

```bash
pip install pytest biopython numpy
```

## Run

From project root:

```bash
pytest tests/ -v
```

## Coverage — 56 Tests

| Module | # | What's Tested |
|:-------|:--|:--------------|
| **`01_parse_and_annotate`** | 6 | DB schema creation, FASTA/CSV loading, HEPN motif detection (`R.{4,6}H`), Protenix JSON generation |
| **`evolution_orchestrator`** | 13 | `compute_fitness` (shift + ipTM + AF2-IG), `EvolutionGym` EMA + bias clipping, mutation extraction (subs/dels/ins), catalytic histidine indices, `build_metadata_override` |
| **`03_pxdesign_wrapper`** | 10 | Frozen REC config, HEPN boundary coords (0-based), `_resolve_unknown_residues` (X→baseline/Gly), `_apply_bias_to_sequence` (linker-only, threshold 0.5) |
| **`hepn_structural_stitch`** | 2 | Wild-type HEPN grafting into designed linker, short-sequence edge case |
| **`pdb_kinematics`** | 7 | Cα HEPN distance (Euclidean), Protenix summary parsing (ipTM/pTM/AF2-IG), missing-file errors, custom chain IDs |
| **`protenix_eval`** | 10 | Mismatch sequence generation (seed reproducibility, exact counts), OFF/ON JSON payloads, off-target JSON, crRNA lookup override |
| **`run_protenix (mocked)`** | 2 | CLI command construction (mini/base tiers), `pred`↔`predict` subcommand fallback |

## Fixtures

Defined in `conftest.py`:

| Fixture | Description |
|:--------|:------------|
| `tmpdir` | Temporary directory (auto-cleanup) |
| `sample_fasta` / `sample_csv` | FASTA + metadata with embedded HEPN motifs |
| `sample_metadata_json` | Variant domain metadata (HEPN coords) |
| `minimal_pdb` | PDB with residues 10 and 50, 30 Å apart |
| `protenix_summary_json` | Mock Protenix summary output |
| `variant_fasta` / `baseline_json` | For mutation extraction and evaluation tests |

## Adding Tests

Follow the existing pattern:
1. Create fixtures in `conftest.py` for reusable test data
2. Use `tmp_path` or `tmpdir` for file I/O
3. Mock external tools (Protenix, PXDesign) — tests must run without GPU
4. Name test files `test_<module>.py` and classes `Test<Feature>`
