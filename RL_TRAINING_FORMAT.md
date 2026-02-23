# RL Training Data Format

The evolution orchestrator appends each variant evaluation to `outputs/rl_gym_data/rl_training_dataset.jsonl`. This data is designed for post-training ProteinMPNN (or discrete diffusion models like DRAKES) to bias future Cas13 designs toward better switches.

## Dataset Format

Each line is a JSON object with:

| Field | Type | Description |
|-------|------|-------------|
| `variant_id` | str | Unique variant name (e.g. `Cas13a_positive_control_variant_3`) |
| `generation` | int | Evolution generation number |
| `baseline_id` | str | Parent baseline ID |
| `crrna_lookup_id` | str | crRNA metadata lookup key |
| `sequence` | str | Designed protein sequence (variant) |
| `baseline_sequence` | str | Parent baseline sequence |
| `mutations` | list[str] | Mutations vs baseline (e.g. `["123_A", "145_G"]`) |
| `fitness` | float | Composite fitness (HEPN shift + ipTM + specificity) |
| `off_dist_A` | float | HEPN distance in OFF state (Å) |
| `on_dist_A` | float | HEPN distance in ON state (Å) |
| `iptm` | float | Protenix ipTM score |
| `af2_ig` | float | AF2-IG confidence score |
| `structure_path` | str \| null | Path to CIF/PDB (ON-state ternary). Null if Protenix failed. |
| `offtarget_by_mismatch` | dict | `{1: dist, 2: dist, 3: dist}` HEPN distance per mismatch count |
| `is_elite` | bool | True if passed elite thresholds (ipTM≥0.85, AF2-IG≥0.80, OFF≥25Å, ON≤12Å) |

## Pipeline Reuse After Post-Training

After fine-tuning ProteinMPNN (or a discrete diffusion inverse-folding model) on this data:

1. **Reward**: Use `fitness` as the training signal (or a derived reward from OFF/ON distances).
2. **Input**: `structure_path` → extract protein backbone; the model learns to generate sequences that fold into that structure and achieve high fitness.
3. **Integration**: Replace or bias the ProteinMPNN step in the PXDesign pipeline with the fine-tuned model. The pipeline flow becomes:
   - PXDesign diffusion → backbone structures
   - **Fine-tuned MPNN** → sequences (biased toward Cas13 switch fitness)
   - Protenix → evaluate OFF/ON, rank

## Conversion for DRAKES

[DRAKES](https://github.com/ChenyuWang-Monica/DRAKES) expects (structure, sequence, reward) for discrete diffusion fine-tuning. Convert with:

- **Structure**: Extract protein chain from `structure_path` (CIF/PDB). Use the designed chain (not crRNA/target).
- **Sequence**: `sequence` field.
- **Reward**: `fitness` (or normalize/scale for stability).

## Suggested Minimum Data

Aim for **500+ records** before post-training. Filter to `is_elite == True` for high-quality supervision, or use all records with fitness-weighted loss.
