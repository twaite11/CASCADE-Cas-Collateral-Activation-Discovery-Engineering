/**
 * F-2 fix: a single source of truth for the HEPN / crRNA color palette.
 *
 * Previously every component (Structure3D, CrrnaSpacer, the legend rendered
 * by VariantDetailDrawer, ...) re-declared its own hex codes, which drifted
 * over time -- HEPN1 was orange in 3D and amber in the legend, the crRNA
 * spacer was rose-300 in the spacer view but #f87171 in old screenshots.
 * Importing from here guarantees parity across every surface.
 *
 * Colors are intentionally tied to Tailwind palette names so designers can
 * tweak `tailwind.config.js` in one place if needed.
 */

/** Hex strings for 3Dmol.js (which expects raw CSS colors, not Tailwind). */
export const PALETTE_HEX = {
  hepn1: "#f97316", // orange-500
  hepn2: "#ef4444", // red-500
  rna: "#38bdf8", // sky-400
  protein: "#6b7280", // zinc-500
  structureBg: "#09090b", // zinc-950
  spacer: "#fb7185", // rose-400
} as const;

/** Tailwind text classes for use in HTML / typography contexts. */
export const PALETTE_TW = {
  hepn1: "text-orange-500",
  hepn2: "text-red-500",
  rna: "text-sky-400",
  protein: "text-zinc-500",
} as const;

/** Per-nucleotide text colors for monospaced crRNA rendering. */
export const NT_COLORS_TW: Record<string, string> = {
  A: "text-emerald-300",
  U: "text-amber-300",
  T: "text-amber-300",
  G: "text-sky-300",
  C: "text-rose-300",
};

/** Used by the structure viewer legend and any other surface that needs the
 *  HEPN/RNA "swatch + label" sequence. */
export const HEPN_LEGEND_ITEMS = [
  { id: "rna", color: PALETTE_HEX.rna, label: "RNA/DNA" },
  { id: "hepn1", color: PALETTE_HEX.hepn1, label: "HEPN1" },
  { id: "hepn2", color: PALETTE_HEX.hepn2, label: "HEPN2" },
  { id: "protein", color: PALETTE_HEX.protein, label: "other" },
] as const;
