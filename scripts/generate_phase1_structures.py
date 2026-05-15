#!/usr/bin/env python3
"""Generate Phase 1 screening structures for confirmed Cas13 baselines.

Runs Protenix mini-tier inference on each baseline JSON in ``jsons/``
and writes the output to ``outputs/phase1_screening/{baseline_id}_pred/``.

This is a one-time bootstrap step; the evolution orchestrator expects
these structures to exist before it can start evolving variants.

Usage (on a GPU machine inside the cascade conda env):
    cd /workspace/CASCADE/scripts
    python generate_phase1_structures.py
"""
import json
import logging
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from utils.protenix_eval import run_protenix_inference  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

JSONS_DIR = os.path.join(os.path.dirname(__file__), "..", "jsons")
PHASE1_DIR = os.path.join(os.path.dirname(__file__), "..", "outputs", "phase1_screening")
VALIDATED_IDS = os.path.join(os.path.dirname(__file__), "..", "outputs", "validated_baseline_ids.txt")


def main() -> int:
    os.makedirs(PHASE1_DIR, exist_ok=True)

    baseline_ids: list[str] = []
    if os.path.isfile(VALIDATED_IDS):
        with open(VALIDATED_IDS) as f:
            baseline_ids = [line.strip() for line in f if line.strip() and not line.startswith("#")]
        log.info("Loaded %d baseline IDs from %s", len(baseline_ids), VALIDATED_IDS)
    else:
        log.warning("No validated_baseline_ids.txt found; scanning jsons/ for all baselines")
        for fn in sorted(os.listdir(JSONS_DIR)):
            if fn.endswith(".json") and not fn.startswith("Cas13a_positive"):
                baseline_ids.append(fn.replace(".json", ""))

    successes = 0
    failures = 0

    for bid in baseline_ids:
        json_path = os.path.join(JSONS_DIR, f"{bid}.json")
        if not os.path.isfile(json_path):
            log.warning("JSON not found for %s, skipping", bid)
            failures += 1
            continue

        pred_dir = os.path.join(PHASE1_DIR, f"{bid}_pred")
        if os.path.isdir(pred_dir) and any(
            f.endswith((".cif", ".pdb")) for f in os.listdir(pred_dir)
        ):
            log.info("Phase 1 structure already exists for %s, skipping", bid)
            successes += 1
            continue

        log.info("Running Protenix mini inference for %s ...", bid)
        try:
            struct_path, summary_path = run_protenix_inference(
                json_path, PHASE1_DIR, model_tier="mini"
            )
            log.info("  -> structure: %s", struct_path)
            log.info("  -> summary:   %s", summary_path)
            successes += 1
        except Exception:
            log.exception("Failed to generate Phase 1 structure for %s", bid)
            failures += 1

    log.info("Phase 1 generation complete: %d succeeded, %d failed", successes, failures)
    return 1 if failures > 0 and successes == 0 else 0


if __name__ == "__main__":
    raise SystemExit(main())
