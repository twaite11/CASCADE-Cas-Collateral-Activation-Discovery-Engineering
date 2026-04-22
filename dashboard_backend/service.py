from __future__ import annotations

import csv
import json
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from statistics import mean
from typing import Any
import time

from .config import DashboardConfig
from .storage import VariantCatalogStore


@dataclass
class DashboardService:
    config: DashboardConfig
    catalog_store: VariantCatalogStore
    cache_ttl_seconds: float = 5.0
    _cache: dict[str, tuple[float, Any]] = field(default_factory=dict)

    def _read_json(self, path: Path, fallback: Any) -> Any:
        if not path.exists():
            return fallback
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except Exception:  # noqa: BLE001
            return fallback

    def _cache_get(self, key: str) -> Any | None:
        now = time.time()
        item = self._cache.get(key)
        if not item:
            return None
        ts, value = item
        if now - ts <= self.cache_ttl_seconds:
            return value
        return None

    def _cache_set(self, key: str, value: Any) -> Any:
        self._cache[key] = (time.time(), value)
        return value

    def invalidate(self) -> None:
        """Drop all cached derived state. Called after a run is promoted so
        the next request re-reads the filesystem."""
        self._cache.clear()

    def _safe_rel(self, path: Path) -> str:
        try:
            return str(path.resolve().relative_to(self.config.cascade_root.resolve()))
        except Exception:  # noqa: BLE001
            return str(path)

    def _read_jsonl(self, path: Path) -> list[dict[str, Any]]:
        if not path.exists():
            return []
        out: list[dict[str, Any]] = []
        for line in path.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                obj = json.loads(line)
                if isinstance(obj, dict):
                    out.append(obj)
            except json.JSONDecodeError:
                continue
        return out

    def _read_csv_rows(self, path: Path) -> list[dict[str, Any]]:
        if not path.exists():
            return []
        rows: list[dict[str, Any]] = []
        try:
            with path.open("r", encoding="utf-8", newline="") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    rows.append(dict(row))
        except Exception:  # noqa: BLE001
            return []
        return rows

    def _index_optimized(self) -> dict[str, dict[str, str]]:
        cached = self._cache_get("index_optimized")
        if cached is not None:
            return cached
        root = self.config.optimized_dir
        if not root.exists():
            return self._cache_set("index_optimized", {})

        out: dict[str, dict[str, str]] = {}
        for fasta in root.glob("*_optimal.fasta"):
            variant_id = fasta.name.replace("_optimal.fasta", "")
            entry = out.setdefault(variant_id, {})
            entry["fasta"] = self._safe_rel(fasta)

        for structure in list(root.glob("*_ternary_complex.pdb")) + list(
            root.glob("*_ternary_complex.cif")
        ):
            suffix = "_ternary_complex" + structure.suffix
            variant_id = structure.name.replace(suffix, "")
            entry = out.setdefault(variant_id, {})
            entry["structure"] = self._safe_rel(structure)

        for crrna in root.glob("*_crRNA.fasta"):
            variant_id = crrna.name.replace("_crRNA.fasta", "")
            entry = out.setdefault(variant_id, {})
            entry["crrna"] = self._safe_rel(crrna)
        return self._cache_set("index_optimized", out)

    def _index_eval_outputs(self) -> dict[str, dict[str, str]]:
        """
        Index nested eval outputs and map to variant and state.

        Expected names include suffixes like:
          {variant_id}_OFF, {variant_id}_ON, {variant_id}_offtarget_...
        """
        cached = self._cache_get("index_eval_outputs")
        if cached is not None:
            return cached
        mapping: dict[str, dict[str, str]] = defaultdict(dict)
        for root in [self.config.fast_eval_dir, self.config.hf_eval_dir]:
            if not root.exists():
                continue

            for pred_dir in root.glob("**/*"):
                if not pred_dir.is_dir():
                    continue
                name = pred_dir.name

                structure_candidates = list(pred_dir.glob("**/*.cif")) + list(pred_dir.glob("**/*.pdb"))
                summary_candidates = list(pred_dir.glob("**/*_summary*.json")) + list(
                    pred_dir.glob("**/*_confidence*.json")
                )
                if not structure_candidates and not summary_candidates:
                    continue

                variant_id = name
                state = "unknown"
                if name.endswith("_OFF"):
                    variant_id = name[: -len("_OFF")]
                    state = "off"
                elif name.endswith("_ON"):
                    variant_id = name[: -len("_ON")]
                    state = "on"
                elif "_offtarget_" in name:
                    variant_id = name.split("_offtarget_")[0]
                    state = "offtarget"

                if structure_candidates:
                    mapping[variant_id][f"{state}_structure"] = self._safe_rel(structure_candidates[0])
                if summary_candidates:
                    mapping[variant_id][f"{state}_summary"] = self._safe_rel(summary_candidates[0])
        return self._cache_set("index_eval_outputs", mapping)

    def _load_validation_index(self) -> dict[str, dict[str, Any]]:
        cached = self._cache_get("validation_index")
        if cached is not None:
            return cached
        report_path = self.config.cascade_root / "outputs" / "repeat_validation_report.csv"
        rows = self._read_csv_rows(report_path)
        idx: dict[str, dict[str, Any]] = {}
        for row in rows:
            key = str(row.get("sequence_id", "")).strip()
            if key:
                idx[key] = row
        return self._cache_set("validation_index", idx)

    def _load_validated_ids(self) -> set[str]:
        cached = self._cache_get("validated_ids")
        if cached is not None:
            return cached
        path = self.config.cascade_root / "outputs" / "validated_baseline_ids.txt"
        if not path.exists():
            return self._cache_set("validated_ids", set())
        out = set()
        for line in path.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if line:
                out.add(line)
        return self._cache_set("validated_ids", out)

    def _threshold_pass(self, row: dict[str, Any]) -> bool:
        iptm = float(row.get("iptm", 0.0) or 0.0)
        af2_ig = float(row.get("af2_ig", 0.0) or 0.0)
        on_dist = float(row.get("on_dist_A", 999.0) or 999.0)
        return (
            iptm >= self.config.min_iptm
            and af2_ig >= self.config.min_af2_ig
            and on_dist <= self.config.max_on_distance
        )

    def _attach_metadata(self, records: list[dict[str, Any]]) -> list[dict[str, Any]]:
        domain_meta = self._read_json(self.config.domain_metadata_path, {})
        catalog_rows = {row.sequence_id: row for row in self.catalog_store.fetch_all_variants()}
        optimized_map = self._index_optimized()
        eval_map = self._index_eval_outputs()
        validation_idx = self._load_validation_index()
        validated_ids = self._load_validated_ids()

        enriched: list[dict[str, Any]] = []
        for row in records:
            variant_id = str(row.get("variant_id", ""))
            baseline_id = str(row.get("baseline_id", ""))
            lookup_id = str(row.get("crrna_lookup_id", "") or baseline_id)

            metadata = domain_meta.get(lookup_id, {})
            catalog = catalog_rows.get(lookup_id)
            optimized_artifacts = optimized_map.get(variant_id, {})
            eval_artifacts = eval_map.get(variant_id, {})
            validation = validation_idx.get(lookup_id)

            is_elite = bool(row.get("is_elite", False))
            threshold_pass = self._threshold_pass(row)
            optimized_switch = is_elite or bool(optimized_artifacts) or threshold_pass

            reasons = []
            if is_elite:
                reasons.append("elite")
            if optimized_artifacts:
                reasons.append("optimized_output")
            if threshold_pass:
                reasons.append("threshold_pass")

            off_dist = float(row.get("off_dist_A", 0.0) or 0.0)
            on_dist = float(row.get("on_dist_A", 0.0) or 0.0)
            enriched.append(
                {
                    **row,
                    "hepn_shift_A": off_dist - on_dist,
                    "optimized_switch": optimized_switch,
                    "optimized_reasons": reasons,
                    "domain_metadata": metadata,
                    "validated_baseline": lookup_id in validated_ids,
                    "validation_metadata": validation,
                    "catalog_metadata": None
                    if catalog is None
                    else {
                        "sra_accession": catalog.sra_accession,
                        "score": catalog.score,
                        "status": catalog.status,
                        "reason": catalog.reason,
                        "hepn1_start": catalog.hepn1_start,
                        "hepn1_end": catalog.hepn1_end,
                        "hepn2_start": catalog.hepn2_start,
                        "hepn2_end": catalog.hepn2_end,
                    },
                    "optimized_artifacts": optimized_artifacts,
                    "eval_artifacts": eval_artifacts,
                }
            )
        return enriched

    def load_variants(self) -> list[dict[str, Any]]:
        cached = self._cache_get("variants")
        if cached is not None:
            return cached
        records = self._read_jsonl(self.config.rl_dataset_path)
        enriched = self._attach_metadata(records)
        out = sorted(
            enriched,
            key=lambda r: (int(r.get("generation", 0)), str(r.get("variant_id", ""))),
            reverse=True,
        )
        return self._cache_set("variants", out)

    def _lineage_index(self, rows: list[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
        idx: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for row in rows:
            lineage = str(row.get("baseline_id", "") or "")
            idx[lineage].append(row)
        for lineage in idx:
            idx[lineage] = sorted(
                idx[lineage],
                key=lambda r: int(r.get("generation", 0) or 0),
            )
        return idx

    def get_variant(self, variant_id: str) -> dict[str, Any] | None:
        all_rows = self.load_variants()
        lineage_idx = self._lineage_index(all_rows)
        for row in all_rows:
            if row.get("variant_id") == variant_id:
                lineage = str(row.get("baseline_id", "") or "")
                lineage_rows = lineage_idx.get(lineage, [])
                lineage_best = (
                    max(lineage_rows, key=lambda x: float(x.get("fitness", -1e9) or -1e9))
                    if lineage_rows
                    else None
                )
                return {
                    **row,
                    "lineage_summary": {
                        "lineage_id": lineage,
                        "lineage_size": len(lineage_rows),
                        "best_variant_id": None if lineage_best is None else lineage_best.get("variant_id"),
                        "best_fitness": None if lineage_best is None else lineage_best.get("fitness"),
                    },
                    "lineage_recent": list(reversed(lineage_rows[-10:])),
                }
        return None

    def get_optimized_summary(self, limit: int = 100) -> list[dict[str, Any]]:
        rows = [r for r in self.load_variants() if r.get("optimized_switch")]
        rows = sorted(rows, key=lambda x: float(x.get("fitness", -1e9) or -1e9), reverse=True)
        out: list[dict[str, Any]] = []
        for row in rows[:limit]:
            out.append(
                {
                    "variant_id": row.get("variant_id"),
                    "baseline_id": row.get("baseline_id"),
                    "generation": row.get("generation"),
                    "fitness": row.get("fitness"),
                    "iptm": row.get("iptm"),
                    "af2_ig": row.get("af2_ig"),
                    "on_dist_A": row.get("on_dist_A"),
                    "off_dist_A": row.get("off_dist_A"),
                    "hepn_shift_A": row.get("hepn_shift_A"),
                    "optimized_reasons": row.get("optimized_reasons", []),
                    "optimized_artifacts": row.get("optimized_artifacts", {}),
                    "eval_artifacts": row.get("eval_artifacts", {}),
                    "domain_metadata": row.get("domain_metadata", {}),
                }
            )
        return out

    def get_overview(self) -> dict[str, Any]:
        variants = self.load_variants()
        generations = [int(v.get("generation", 0) or 0) for v in variants]
        fitness_values = [float(v.get("fitness", 0.0) or 0.0) for v in variants]
        iptm_values = [float(v.get("iptm", 0.0) or 0.0) for v in variants]
        af2_values = [float(v.get("af2_ig", 0.0) or 0.0) for v in variants]
        optimized_count = sum(1 for v in variants if v.get("optimized_switch"))
        validated_count = sum(1 for v in variants if v.get("validated_baseline"))

        now_iso = datetime.now(timezone.utc).isoformat()
        dataset_mtime = (
            datetime.fromtimestamp(self.config.rl_dataset_path.stat().st_mtime, tz=timezone.utc).isoformat()
            if self.config.rl_dataset_path.exists()
            else None
        )
        stale_data = False
        if self.config.rl_dataset_path.exists():
            age_seconds = datetime.now(timezone.utc).timestamp() - self.config.rl_dataset_path.stat().st_mtime
            stale_data = age_seconds > 3600 * 6

        warnings: list[str] = []
        if not self.config.rl_dataset_path.exists():
            warnings.append("rl_dataset_missing")
        if not self.config.sqlite_db_path.exists():
            warnings.append("sqlite_catalog_missing")
        if stale_data:
            warnings.append("rl_dataset_stale_over_6h")
        required_fields = {
            "variant_id",
            "generation",
            "baseline_id",
            "fitness",
            "off_dist_A",
            "on_dist_A",
            "iptm",
            "af2_ig",
        }
        if variants:
            missing_count = sum(
                1
                for row in variants
                if any(k not in row or row.get(k) is None for k in required_fields)
            )
            if missing_count > 0:
                warnings.append(f"schema_drift_missing_required_fields:{missing_count}")

        return {
            "generated_at": now_iso,
            "dataset_last_updated": dataset_mtime,
            "total_variants": len(variants),
            "max_generation": max(generations) if generations else 0,
            "optimized_switches": optimized_count,
            "validated_baseline_records": validated_count,
            "elite_count": sum(1 for v in variants if v.get("is_elite")),
            "mean_fitness": mean(fitness_values) if fitness_values else 0.0,
            "mean_iptm": mean(iptm_values) if iptm_values else 0.0,
            "mean_af2_ig": mean(af2_values) if af2_values else 0.0,
            "generation_histogram": dict(Counter(generations)),
            "warnings": warnings,
        }

    def get_pipeline_health(self) -> dict[str, Any]:
        variants = self.load_variants()
        if not variants:
            return {
                "total_records": 0,
                "structure_success_rate": 0.0,
                "full_ternary_rate": 0.0,
                "failure_count": 0,
                "best_fitness": None,
                "stage_stats": {},
                "throughput_variants_per_hour": 0.0,
            }

        structure_ok = sum(1 for v in variants if v.get("structure_path"))
        full_ternary = sum(1 for v in variants if float(v.get("iptm", 0.0) or 0.0) >= self.config.min_iptm)
        failures = sum(1 for v in variants if float(v.get("on_dist_A", 999.0) or 999.0) >= 900.0)
        best = max(variants, key=lambda x: float(x.get("fitness", -1e9) or -1e9))

        stage_stats = {
            "pxdesign_fail_like": failures,
            "mini_eval_success": sum(
                1
                for v in variants
                if float(v.get("off_dist_A", 999.0) or 999.0) < 900
                and float(v.get("on_dist_A", 999.0) or 999.0) < 900
            ),
            "base_eval_success": sum(
                1 for v in variants if float(v.get("iptm", 0.0) or 0.0) >= self.config.min_iptm
            ),
            "msa_reuse_observed": self._count_log_keyword(r"Reusing cached MSA output"),
        }

        throughput_variants_per_hour = 0.0
        if self.config.rl_dataset_path.exists() and len(variants) > 0:
            first = self.config.rl_dataset_path.stat().st_mtime
            logs = sorted((self.config.cascade_root / "logs").glob("evolution_*.log"))
            if logs:
                first = min(first, logs[0].stat().st_mtime)
            elapsed_h = max((datetime.now(timezone.utc).timestamp() - first) / 3600.0, 1e-6)
            throughput_variants_per_hour = len(variants) / elapsed_h

        return {
            "total_records": len(variants),
            "structure_success_rate": structure_ok / len(variants),
            "full_ternary_rate": full_ternary / len(variants),
            "failure_count": failures,
            "best_variant_id": best.get("variant_id"),
            "best_fitness": float(best.get("fitness", 0.0) or 0.0),
            "best_iptm": float(best.get("iptm", 0.0) or 0.0),
            "best_af2_ig": float(best.get("af2_ig", 0.0) or 0.0),
            "stage_stats": stage_stats,
            "throughput_variants_per_hour": throughput_variants_per_hour,
        }

    def _count_log_keyword(self, pattern: str) -> int:
        logs_dir = self.config.cascade_root / "logs"
        if not logs_dir.exists():
            return 0
        rx = re.compile(pattern)
        count = 0
        for log_file in sorted(logs_dir.glob("evolution_*.log"))[-3:]:
            try:
                text = log_file.read_text(encoding="utf-8", errors="ignore")
            except Exception:  # noqa: BLE001
                continue
            count += len(rx.findall(text))
        return count

    def ping(self) -> dict[str, Any]:
        db = self.catalog_store.ping()
        return {
            "status": "ok",
            "db": db,
            "paths": {
                "cascade_root": str(self.config.cascade_root),
                "rl_dataset_path": str(self.config.rl_dataset_path),
            },
        }

    def get_production_summary(self) -> dict[str, Any]:
        rows = self.load_variants()
        by_lineage: dict[str, dict[str, Any]] = defaultdict(
            lambda: {"generated": 0, "optimized": 0, "elite": 0}
        )
        for r in rows:
            lineage = str(r.get("baseline_id", "") or "")
            by_lineage[lineage]["generated"] += 1
            if r.get("optimized_switch"):
                by_lineage[lineage]["optimized"] += 1
            if r.get("is_elite"):
                by_lineage[lineage]["elite"] += 1

        lineages = []
        for lineage, counts in by_lineage.items():
            generated = counts["generated"] or 1
            lineages.append(
                {
                    "lineage_id": lineage,
                    **counts,
                    "optimized_yield": counts["optimized"] / generated,
                    "elite_yield": counts["elite"] / generated,
                }
            )
        lineages.sort(key=lambda x: x["optimized_yield"], reverse=True)

        total = len(rows) or 1
        passed_filter = sum(1 for r in rows if float(r.get("on_dist_A", 999.0) or 999.0) <= self.config.max_on_distance)
        optimized = sum(1 for r in rows if r.get("optimized_switch"))
        elite = sum(1 for r in rows if r.get("is_elite"))
        return {
            "totals": {
                "generated": len(rows),
                "passed_filter": passed_filter,
                "optimized": optimized,
                "elite": elite,
            },
            "funnel": {
                "generated_to_passed": passed_filter / total,
                "passed_to_optimized": optimized / max(passed_filter, 1),
                "optimized_to_elite": elite / max(optimized, 1),
            },
            "lineages": lineages[:50],
        }
