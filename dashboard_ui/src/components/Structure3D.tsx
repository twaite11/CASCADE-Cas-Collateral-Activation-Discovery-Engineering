import { useEffect, useRef, useState } from "react";
import * as $3Dmol from "3dmol";

import { structureFileUrl, withAuth } from "@/lib/api";
import type { DomainMetadata } from "@/lib/types";
import { HEPN_LEGEND_ITEMS, PALETTE_HEX } from "@/lib/palette";

interface Props {
  path: string;
  hepn?: DomainMetadata;
  heightPx?: number;
}

/**
 * 3Dmol.js cartoon viewer. Loads a CIF/PDB via /api/structure-file, colors
 * HEPN1/HEPN2 ranges if domain metadata is provided, tints everything else
 * neutral grey, and renders RNA/DNA chains as cartoon ribbons.
 */
export function Structure3D({ path, hepn, heightPx = 280 }: Props) {
  const hostRef = useRef<HTMLDivElement | null>(null);
  const viewerRef = useRef<any>(null);
  const [error, setError] = useState<string | null>(null);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    let disposed = false;
    // Phase-3 fix: abortable fetch so an unmount mid-load doesn't leak a
    // half-completed render into the next mount.
    const ctrl = new AbortController();

    async function boot() {
      if (!hostRef.current) return;

      hostRef.current.innerHTML = "";
      setLoading(true);
      setError(null);

      let data: string;
      let fmt: "cif" | "pdb";
      try {
        const resp = await fetch(structureFileUrl(path), withAuth({ signal: ctrl.signal }));
        if (!resp.ok) throw new Error(`${resp.status} ${resp.statusText}`);
        data = await resp.text();
        fmt = path.toLowerCase().endsWith(".cif") ? "cif" : "pdb";
      } catch (e: any) {
        if (!disposed) {
          // AbortError is expected on rapid prop changes; don't surface it.
          if (e?.name !== "AbortError") {
            setError(String(e?.message ?? e));
            setLoading(false);
          }
        }
        return;
      }

      if (disposed || !hostRef.current) return;

      const viewer = $3Dmol.createViewer(hostRef.current, {
        backgroundColor: PALETTE_HEX.structureBg,
      });
      viewerRef.current = viewer;

      viewer.addModel(data, fmt);

      // Neutral base cartoon for everything
      viewer.setStyle(
        {},
        { cartoon: { color: PALETTE_HEX.protein, thickness: 0.6 } },
      );
      // RNA / DNA: blue
      viewer.setStyle(
        { resn: ["A", "G", "C", "U", "T", "DA", "DG", "DC", "DT"] },
        { cartoon: { color: PALETTE_HEX.rna, thickness: 0.8 } },
      );

      // HEPN coloring — palette is centralised in @/lib/palette so the 3D
      // viewer, the spacer track, and the legend never drift apart (F-2).
      if (hepn?.hepn1_start != null && hepn?.hepn1_end != null) {
        viewer.addStyle(
          { resi: `${hepn.hepn1_start}-${hepn.hepn1_end}` },
          { cartoon: { color: PALETTE_HEX.hepn1, thickness: 1.0 } },
        );
      }
      if (hepn?.hepn2_start != null && hepn?.hepn2_end != null) {
        viewer.addStyle(
          { resi: `${hepn.hepn2_start}-${hepn.hepn2_end}` },
          { cartoon: { color: PALETTE_HEX.hepn2, thickness: 1.0 } },
        );
      }

      viewer.zoomTo();
      viewer.render();
      viewer.zoom(0.9, 500);

      if (!disposed) setLoading(false);
    }

    boot();

    return () => {
      disposed = true;
      ctrl.abort();
      try {
        viewerRef.current?.clear();
        viewerRef.current?.removeAllModels();
      } catch {
        /* noop */
      }
      viewerRef.current = null;
    };
  }, [path, hepn?.hepn1_start, hepn?.hepn1_end, hepn?.hepn2_start, hepn?.hepn2_end]);

  return (
    <div className="relative">
      <div
        ref={hostRef}
        className="overflow-hidden rounded-md border bg-zinc-950"
        style={{ height: heightPx, position: "relative" }}
      />
      {loading && (
        <div className="pointer-events-none absolute inset-0 flex items-center justify-center text-xs text-muted-foreground">
          loading structure…
        </div>
      )}
      {error && (
        <div className="pointer-events-none absolute inset-0 flex items-center justify-center text-xs text-red-400">
          {error}
        </div>
      )}
      <Legend hepn={hepn} />
    </div>
  );
}

function Legend({ hepn }: { hepn?: DomainMetadata }) {
  // F-2 fix: legend reads from the shared palette + filters out HEPN dots
  // that aren't actually rendered for this variant.
  const hasHepn1 = hepn?.hepn1_start != null;
  const hasHepn2 = hepn?.hepn2_start != null;
  return (
    <div className="mt-1 flex flex-wrap gap-x-3 gap-y-1 text-[10px] text-muted-foreground">
      {HEPN_LEGEND_ITEMS.filter((it) => {
        if (it.id === "hepn1") return hasHepn1;
        if (it.id === "hepn2") return hasHepn2;
        return true;
      }).map((it) => (
        <LegendDot key={it.id} color={it.color} label={it.label} />
      ))}
    </div>
  );
}

function LegendDot({ color, label }: { color: string; label: string }) {
  return (
    <span className="inline-flex items-center gap-1">
      <span
        className="inline-block h-2 w-2 rounded-full"
        style={{ background: color }}
      />
      {label}
    </span>
  );
}
