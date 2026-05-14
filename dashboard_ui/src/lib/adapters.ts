/**
 * C-10/C-11 fix: a single, type-safe normalizer for the variant blobs the
 * backend returns from `/api/variant/{id}` and `/api/optimized-switches`.
 *
 * The backend mixes three data sources together inside one record:
 *   - `domain_metadata`  (variant_domain_metadata.json, by `crrna_lookup_id`)
 *   - `catalog_metadata` (cas13_variants.db row, by `crrna_lookup_id`)
 *   - top-level fields  (RL dataset jsonl)
 *
 * Before this adapter every component reached into those blobs with
 * `(v.domain_metadata as any)?.crrna_repeat ?? (v.catalog_metadata as any)?...`
 * pattern, which silently fell through to "null" when one source moved keys.
 *
 * This module pins the contract: read all three sources in a known order,
 * coerce to a strongly-typed `NormalizedVariant`, and let the caller use
 * `nv.hepn`, `nv.crrna.repeat`, etc. with no defensive casts.
 */
import type { DomainMetadata, OptimizedVariant } from "./types";

export interface NormalizedHepn {
  hepn1Start: number | null;
  hepn1End: number | null;
  hepn2Start: number | null;
  hepn2End: number | null;
  /** True iff at least one HEPN range is fully resolved. */
  hasRange: boolean;
}

export interface NormalizedCrrna {
  repeat: string | null;
  spacer: string | null;
  lookupId: string | null;
}

export interface NormalizedArtifacts {
  optimized: {
    fasta: string | null;
    structure: string | null;
    crrna: string | null;
  };
  evals: {
    on: { structure: string | null; summary: string | null };
    off: { structure: string | null; summary: string | null };
    offtarget: { structure: string | null; summary: string | null };
  };
  /** Best-effort: the most representative structure to show by default. */
  primaryStructure: string | null;
}

export interface NormalizedVariant {
  variantId: string;
  baselineId: string | null;
  generation: number | null;
  fitness: number | null;
  iptm: number | null;
  af2Ig: number | null;
  onDistA: number | null;
  offDistA: number | null;
  hepnShiftA: number | null;
  hepn: NormalizedHepn;
  crrna: NormalizedCrrna;
  artifacts: NormalizedArtifacts;
  optimizedReasons: string[];
  validatedBaseline: boolean;
  isElite: boolean;
  /** Untouched original blob, for the "raw JSON" drawer section. */
  raw: Record<string, unknown>;
}

function asNumberOrNull(v: unknown): number | null {
  if (v === null || v === undefined || v === "") return null;
  const n = Number(v);
  return Number.isFinite(n) ? n : null;
}

function asStringOrNull(v: unknown): string | null {
  if (v === null || v === undefined) return null;
  if (typeof v === "string") return v.trim() || null;
  return String(v).trim() || null;
}

function asStringArray(v: unknown): string[] {
  if (!Array.isArray(v)) return [];
  return v.filter((x): x is string => typeof x === "string" && x.length > 0);
}

/**
 * Pick the first non-null value from a list of source objects under a list
 * of possible keys.  This is the "alias resolver" that consolidates the
 * many ways the backend can spell the same field across the three blobs.
 */
function pickStr(
  sources: Array<Record<string, unknown> | undefined | null>,
  keys: string[],
): string | null {
  for (const src of sources) {
    if (!src) continue;
    for (const k of keys) {
      const v = src[k];
      const out = asStringOrNull(v);
      if (out !== null) return out;
    }
  }
  return null;
}

function pickNum(
  sources: Array<Record<string, unknown> | undefined | null>,
  keys: string[],
): number | null {
  for (const src of sources) {
    if (!src) continue;
    for (const k of keys) {
      const v = src[k];
      const out = asNumberOrNull(v);
      if (out !== null) return out;
    }
  }
  return null;
}

/**
 * Normalize a variant-detail (or optimized-switches row) payload.  The
 * input is intentionally typed as `unknown`/`Record<string, unknown>` so the
 * caller doesn't accidentally assume the strict OptimizedVariant shape.
 */
export function normalizeVariant(v: Record<string, unknown>): NormalizedVariant {
  const domain = (v.domain_metadata as Record<string, unknown>) ?? {};
  const catalog = (v.catalog_metadata as Record<string, unknown>) ?? {};
  const optArtifactsBlob =
    (v.optimized_artifacts as Record<string, unknown>) ?? {};
  const evalArtifactsBlob =
    (v.eval_artifacts as Record<string, unknown>) ?? {};

  // HEPN: prefer domain_metadata, fall back to catalog_metadata.
  const hepn: NormalizedHepn = {
    hepn1Start: pickNum([domain, catalog], ["hepn1_start"]),
    hepn1End: pickNum([domain, catalog], ["hepn1_end"]),
    hepn2Start: pickNum([domain, catalog], ["hepn2_start"]),
    hepn2End: pickNum([domain, catalog], ["hepn2_end"]),
    hasRange: false,
  };
  hepn.hasRange =
    (hepn.hepn1Start !== null && hepn.hepn1End !== null) ||
    (hepn.hepn2Start !== null && hepn.hepn2End !== null);

  // crRNA: backend has used both `crrna_repeat`/`crrna_spacer` (Baseline,
  // catalog metadata) and `crRNA_repeat`/`crRNA_spacer` (older
  // domain_metadata exports).  Accept both.
  const crrna: NormalizedCrrna = {
    repeat: pickStr(
      [domain, catalog, v],
      ["crrna_repeat", "crRNA_repeat", "dr", "direct_repeat", "repeat"],
    ),
    spacer: pickStr(
      [domain, catalog, v],
      ["crrna_spacer", "crRNA_spacer", "spacer"],
    ),
    lookupId: pickStr([v, catalog], ["crrna_lookup_id", "baseline_id"]),
  };

  // Artifacts: the structure of these blobs is well-defined server-side
  // but we coerce to strings + null defensively so a typo in the eval
  // pipeline can't crash <Structure3D/>.
  const optimized = {
    fasta: asStringOrNull(optArtifactsBlob.fasta),
    structure: asStringOrNull(optArtifactsBlob.structure),
    crrna: asStringOrNull(optArtifactsBlob.crrna),
  };
  const evals = {
    on: {
      structure: asStringOrNull(evalArtifactsBlob.on_structure),
      summary: asStringOrNull(evalArtifactsBlob.on_summary),
    },
    off: {
      structure: asStringOrNull(evalArtifactsBlob.off_structure),
      summary: asStringOrNull(evalArtifactsBlob.off_summary),
    },
    offtarget: {
      structure: asStringOrNull(evalArtifactsBlob.offtarget_structure),
      summary: asStringOrNull(evalArtifactsBlob.offtarget_summary),
    },
  };
  const primaryStructure =
    optimized.structure ?? evals.on.structure ?? evals.off.structure ?? null;

  return {
    variantId: asStringOrNull(v.variant_id) ?? "",
    baselineId: asStringOrNull(v.baseline_id),
    generation: asNumberOrNull(v.generation),
    fitness: asNumberOrNull(v.fitness),
    iptm: asNumberOrNull(v.iptm),
    af2Ig: asNumberOrNull(v.af2_ig),
    onDistA: asNumberOrNull(v.on_dist_A),
    offDistA: asNumberOrNull(v.off_dist_A),
    hepnShiftA: asNumberOrNull(v.hepn_shift_A),
    hepn,
    crrna,
    artifacts: { optimized, evals, primaryStructure },
    optimizedReasons: asStringArray(v.optimized_reasons),
    validatedBaseline: Boolean(v.validated_baseline),
    isElite: Boolean(v.is_elite),
    raw: v,
  };
}

/** Subset projection used by the structure viewer. */
export function toDomainMetadata(hepn: NormalizedHepn): DomainMetadata {
  return {
    hepn1_start: hepn.hepn1Start ?? undefined,
    hepn1_end: hepn.hepn1End ?? undefined,
    hepn2_start: hepn.hepn2Start ?? undefined,
    hepn2_end: hepn.hepn2End ?? undefined,
  };
}

/** Adapter so OptimizedVariant rows from the sidebar can flow through the
 *  same normalizer pipeline as variant-detail blobs. */
export function normalizeOptimizedRow(v: OptimizedVariant): NormalizedVariant {
  return normalizeVariant(v as unknown as Record<string, unknown>);
}
