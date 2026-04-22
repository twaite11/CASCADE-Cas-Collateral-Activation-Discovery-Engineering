import { useEffect, useMemo } from "react";
import { useQuery } from "@tanstack/react-query";
import {
  X,
  Download,
  Copy,
  GitBranch,
  Dna,
  Box,
  FileText,
} from "lucide-react";

import { api, structureFileUrl } from "@/lib/api";
import type { DomainMetadata } from "@/lib/types";
import { cn } from "@/lib/utils";
import { toast } from "@/components/ui/toast";
import { useVariantDrawer } from "@/stores/variantDrawer";

import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Structure3D } from "@/components/Structure3D";
import { CrrnaSpacer } from "@/components/CrrnaSpacer";

export function VariantDetailDrawer() {
  const { variantId, compareWith, close, setCompareWith } = useVariantDrawer();
  const open = variantId !== null;

  useEffect(() => {
    if (!open) return;
    const onKey = (e: KeyboardEvent) => {
      if (e.key === "Escape") close();
    };
    window.addEventListener("keydown", onKey);
    document.body.style.overflow = "hidden";
    return () => {
      window.removeEventListener("keydown", onKey);
      document.body.style.overflow = "";
    };
  }, [open, close]);

  if (!open) return null;

  return (
    <div className="fixed inset-0 z-50 flex">
      <div
        className="flex-1 bg-black/60 backdrop-blur-sm"
        onClick={close}
        aria-hidden
      />
      <aside className="flex h-full w-full max-w-3xl flex-col border-l bg-background shadow-2xl">
        <DrawerContent
          variantId={variantId!}
          compareWith={compareWith}
          onClose={close}
          onCompareChange={setCompareWith}
        />
      </aside>
    </div>
  );
}

function DrawerContent({
  variantId,
  compareWith,
  onClose,
  onCompareChange,
}: {
  variantId: string;
  compareWith: string | null;
  onClose: () => void;
  onCompareChange: (id: string | null) => void;
}) {
  const detailQ = useQuery({
    queryKey: ["variant", variantId],
    queryFn: () => api.variantDetail(variantId),
    staleTime: 60_000,
  });

  const v = detailQ.data ?? {};
  const hepn = (v.domain_metadata ?? {}) as DomainMetadata;
  const optArtifacts = (v.optimized_artifacts ?? {}) as Record<string, string>;
  const evalArtifacts = (v.eval_artifacts ?? {}) as Record<string, string>;
  const primaryStructure =
    optArtifacts.structure ??
    evalArtifacts.on_structure ??
    evalArtifacts.off_structure ??
    null;

  const lineage = useMemo(() => extractLineage(v), [v]);

  function copyId() {
    navigator.clipboard.writeText(variantId).catch(() => undefined);
    toast({ kind: "info", title: "Variant ID copied" });
  }

  return (
    <>
      <header className="flex items-center justify-between border-b px-6 py-4">
        <div className="min-w-0">
          <div className="flex items-center gap-2">
            <h2 className="truncate font-mono text-lg font-semibold">
              {variantId}
            </h2>
            <Button size="icon" variant="ghost" onClick={copyId} title="Copy ID">
              <Copy className="h-3.5 w-3.5" />
            </Button>
          </div>
          <p className="truncate text-xs text-muted-foreground">
            {String(v.baseline_id ?? "—")} · gen {String(v.generation ?? "—")}
          </p>
        </div>
        <Button variant="ghost" size="icon" onClick={onClose} title="Close (Esc)">
          <X className="h-5 w-5" />
        </Button>
      </header>

      <div className="flex-1 overflow-y-auto px-6 py-4">
        {detailQ.isLoading && (
          <p className="text-sm text-muted-foreground">loading variant…</p>
        )}
        {detailQ.error && (
          <p className="text-sm text-red-400">
            failed: {(detailQ.error as Error).message}
          </p>
        )}
        {detailQ.isSuccess && (
          <div className="space-y-5">
            <MetricsGrid v={v} />

            <Card>
              <CardHeader className="pb-2">
                <CardTitle className="flex items-center gap-2 text-sm">
                  <GitBranch className="h-4 w-4" /> Lineage
                </CardTitle>
              </CardHeader>
              <CardContent>
                {lineage.length === 0 ? (
                  <p className="text-xs text-muted-foreground">
                    No parent recorded (likely a root baseline).
                  </p>
                ) : (
                  <ol className="space-y-1 font-mono text-xs">
                    {lineage.map((id, i) => (
                      <li key={id} className="flex items-center gap-2">
                        <span
                          className={cn(
                            "inline-flex h-5 w-5 items-center justify-center rounded-full border text-[10px]",
                            i === lineage.length - 1
                              ? "border-primary text-primary"
                              : "border-muted-foreground/30 text-muted-foreground",
                          )}
                        >
                          {i}
                        </span>
                        {id}
                      </li>
                    ))}
                  </ol>
                )}
              </CardContent>
            </Card>

            <Card>
              <CardHeader className="pb-2">
                <CardTitle className="flex items-center gap-2 text-sm">
                  <Dna className="h-4 w-4" /> crRNA
                </CardTitle>
              </CardHeader>
              <CardContent>
                <CrrnaSpacer
                  repeat={
                    (v.domain_metadata as any)?.crrna_repeat ??
                    (v.catalog_metadata as any)?.crrna_repeat ??
                    null
                  }
                  spacer={
                    (v.domain_metadata as any)?.crrna_spacer ??
                    (v.catalog_metadata as any)?.crrna_spacer ??
                    null
                  }
                />
              </CardContent>
            </Card>

            {primaryStructure && (
              <Card>
                <CardHeader className="flex-row items-center justify-between space-y-0 pb-2">
                  <CardTitle className="flex items-center gap-2 text-sm">
                    <Box className="h-4 w-4" /> 3D structure
                  </CardTitle>
                  <div className="flex items-center gap-1">
                    <label className="text-[11px] text-muted-foreground">
                      Compare with:
                    </label>
                    <input
                      type="text"
                      value={compareWith ?? ""}
                      placeholder="other variant id"
                      onChange={(e) =>
                        onCompareChange(e.target.value.trim() || null)
                      }
                      className="h-7 w-44 rounded-md border bg-background px-2 text-xs font-mono"
                    />
                  </div>
                </CardHeader>
                <CardContent>
                  {compareWith ? (
                    <div className="grid grid-cols-2 gap-3">
                      <div>
                        <p className="mb-1 text-[10px] text-muted-foreground">
                          {variantId}
                        </p>
                        <Structure3D
                          path={primaryStructure}
                          hepn={hepn}
                          heightPx={260}
                        />
                      </div>
                      <CompareStructure id={compareWith} />
                    </div>
                  ) : (
                    <Structure3D
                      path={primaryStructure}
                      hepn={hepn}
                      heightPx={320}
                    />
                  )}
                </CardContent>
              </Card>
            )}

            <ArtifactsCard
              optArtifacts={optArtifacts}
              evalArtifacts={evalArtifacts}
            />

            <RawMetadataCard data={v} />
          </div>
        )}
      </div>
    </>
  );
}

function CompareStructure({ id }: { id: string }) {
  const q = useQuery({
    queryKey: ["variant", id],
    queryFn: () => api.variantDetail(id),
    staleTime: 60_000,
  });
  const v = q.data ?? {};
  const hepn = (v.domain_metadata ?? {}) as DomainMetadata;
  const optArtifacts = (v.optimized_artifacts ?? {}) as Record<string, string>;
  const evalArtifacts = (v.eval_artifacts ?? {}) as Record<string, string>;
  const path =
    optArtifacts.structure ??
    evalArtifacts.on_structure ??
    evalArtifacts.off_structure ??
    null;

  return (
    <div>
      <p className="mb-1 text-[10px] text-muted-foreground">{id}</p>
      {q.isLoading && (
        <div className="flex h-64 items-center justify-center rounded-md border bg-zinc-950 text-xs text-muted-foreground">
          loading…
        </div>
      )}
      {q.isError && (
        <div className="flex h-64 items-center justify-center rounded-md border bg-zinc-950 text-xs text-red-400">
          {(q.error as Error).message}
        </div>
      )}
      {q.isSuccess && path && (
        <Structure3D path={path} hepn={hepn} heightPx={260} />
      )}
      {q.isSuccess && !path && (
        <div className="flex h-64 items-center justify-center rounded-md border bg-zinc-950 text-xs text-muted-foreground">
          no structure for {id}
        </div>
      )}
    </div>
  );
}

function MetricsGrid({ v }: { v: Record<string, unknown> }) {
  const cells: Array<{ label: string; value: unknown; unit?: string }> = [
    { label: "fitness", value: v.fitness },
    { label: "ipTM", value: v.iptm },
    { label: "AF2-IG", value: v.af2_ig },
    { label: "ON dist", value: v.on_dist_A, unit: "Å" },
    { label: "OFF dist", value: v.off_dist_A, unit: "Å" },
    { label: "HEPN Δ", value: v.hepn_shift_A, unit: "Å" },
  ];
  return (
    <Card>
      <CardContent className="grid grid-cols-6 gap-2 p-3 text-center">
        {cells.map((c) => (
          <div key={c.label} className="rounded bg-muted/40 p-2">
            <div className="text-[10px] uppercase text-muted-foreground">
              {c.label}
            </div>
            <div className="font-mono text-sm">
              {formatNum(c.value)}
              {c.unit && (
                <span className="ml-0.5 text-[10px] text-muted-foreground">
                  {c.unit}
                </span>
              )}
            </div>
          </div>
        ))}
      </CardContent>
    </Card>
  );
}

function ArtifactsCard({
  optArtifacts,
  evalArtifacts,
}: {
  optArtifacts: Record<string, string>;
  evalArtifacts: Record<string, string>;
}) {
  const links: Array<{ label: string; path: string }> = [];
  if (optArtifacts.fasta)
    links.push({ label: "optimized FASTA", path: optArtifacts.fasta });
  if (optArtifacts.structure)
    links.push({ label: "optimized structure", path: optArtifacts.structure });
  if (optArtifacts.crrna)
    links.push({ label: "crRNA FASTA", path: optArtifacts.crrna });
  for (const [k, path] of Object.entries(evalArtifacts)) {
    links.push({ label: k.replace(/_/g, " "), path });
  }

  if (links.length === 0) return null;

  return (
    <Card>
      <CardHeader className="pb-2">
        <CardTitle className="flex items-center gap-2 text-sm">
          <Download className="h-4 w-4" /> Artifacts
        </CardTitle>
      </CardHeader>
      <CardContent className="grid grid-cols-2 gap-1">
        {links.map(({ label, path }) => (
          <a
            key={label + path}
            href={structureFileUrl(path)}
            download
            className="flex items-center justify-between gap-2 rounded border bg-muted/20 px-2 py-1.5 text-xs hover:bg-muted/40"
            title={path}
          >
            <span className="truncate">{label}</span>
            <Download className="h-3 w-3 shrink-0 text-muted-foreground" />
          </a>
        ))}
      </CardContent>
    </Card>
  );
}

function RawMetadataCard({ data }: { data: Record<string, unknown> }) {
  return (
    <Card>
      <CardHeader className="pb-2">
        <CardTitle className="flex items-center gap-2 text-sm">
          <FileText className="h-4 w-4" /> Raw metadata
        </CardTitle>
      </CardHeader>
      <CardContent>
        <div className="flex flex-wrap gap-1">
          {Array.isArray((data as any).optimized_reasons) &&
            ((data as any).optimized_reasons as string[]).map((r) => (
              <Badge key={r} variant="outline">
                {r}
              </Badge>
            ))}
        </div>
        <details className="mt-2">
          <summary className="cursor-pointer text-xs text-muted-foreground hover:text-foreground">
            full JSON
          </summary>
          <pre className="mt-2 max-h-80 overflow-auto rounded border bg-zinc-950 p-2 text-[10px] leading-snug text-zinc-300">
            {JSON.stringify(data, null, 2)}
          </pre>
        </details>
      </CardContent>
    </Card>
  );
}

function extractLineage(v: Record<string, unknown>): string[] {
  const parents = (v.parent_ids ?? v.lineage ?? v.ancestry) as
    | string[]
    | undefined;
  if (Array.isArray(parents)) return parents.filter(Boolean);
  const parent = (v.parent_id ?? v.parent) as string | undefined;
  if (parent) return [parent];
  return [];
}

function formatNum(n: unknown): string {
  if (n == null || n === "") return "—";
  const num = Number(n);
  if (Number.isNaN(num)) return String(n);
  if (Math.abs(num) >= 100) return num.toFixed(1);
  return num.toFixed(3);
}
