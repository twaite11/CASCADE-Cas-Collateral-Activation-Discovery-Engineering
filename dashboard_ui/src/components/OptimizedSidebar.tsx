import { useMemo, useState } from "react";
import { useQuery } from "@tanstack/react-query";
import {
  Sparkles,
  PanelRightClose,
  PanelRightOpen,
  Filter,
  ChevronDown,
  ChevronRight,
  Box,
} from "lucide-react";

import { api } from "@/lib/api";
import type { OptimizedVariant } from "@/lib/types";
import { cn } from "@/lib/utils";

import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Structure3D } from "@/components/Structure3D";
import { useVariantDrawer } from "@/stores/variantDrawer";

type StructureState = "on" | "off" | "offtarget" | "optimized";

export function OptimizedSidebar() {
  const [collapsed, setCollapsed] = useState(false);
  const [minFitness, setMinFitness] = useState<number | "">("");
  const [minIptm, setMinIptm] = useState<number | "">("");
  const [generation, setGeneration] = useState<number | "">("");
  const [lineage, setLineage] = useState("");

  const { data, isLoading, error } = useQuery({
    queryKey: ["optimized-switches"],
    queryFn: () => api.optimized(500),
    refetchInterval: 15_000,
  });

  const rows = data?.rows ?? [];

  const lineages = useMemo(() => {
    const set = new Set<string>();
    for (const r of rows) if (r.baseline_id) set.add(r.baseline_id);
    return Array.from(set).sort();
  }, [rows]);

  const filtered = useMemo(() => {
    return rows.filter((r) => {
      if (minFitness !== "" && (r.fitness ?? -Infinity) < minFitness) return false;
      if (minIptm !== "" && (r.iptm ?? 0) < minIptm) return false;
      if (generation !== "" && r.generation !== generation) return false;
      if (lineage && r.baseline_id !== lineage) return false;
      return true;
    });
  }, [rows, minFitness, minIptm, generation, lineage]);

  if (collapsed) {
    return (
      <aside className="hidden w-10 shrink-0 border-l bg-card/30 p-2 xl:flex xl:flex-col xl:items-center">
        <Button
          variant="ghost"
          size="icon"
          onClick={() => setCollapsed(false)}
          title="Expand Optimized Switches"
        >
          <PanelRightOpen className="h-4 w-4" />
        </Button>
        <Sparkles className="mt-4 h-4 w-4 text-muted-foreground" />
      </aside>
    );
  }

  return (
    <aside className="hidden w-[420px] shrink-0 flex-col overflow-hidden border-l bg-card/30 xl:flex">
      <header className="flex items-center justify-between border-b px-4 py-3">
        <div className="flex items-center gap-2">
          <Sparkles className="h-4 w-4 text-amber-400" />
          <h2 className="text-sm font-semibold">Optimized Switches</h2>
          <Badge variant="outline" className="text-[10px]">
            {filtered.length}/{rows.length}
          </Badge>
        </div>
        <Button
          variant="ghost"
          size="icon"
          onClick={() => setCollapsed(true)}
          title="Collapse"
        >
          <PanelRightClose className="h-4 w-4" />
        </Button>
      </header>

      <div className="border-b bg-muted/20 px-4 py-3">
        <div className="mb-2 flex items-center gap-1 text-[11px] font-medium uppercase text-muted-foreground">
          <Filter className="h-3 w-3" /> Filters
        </div>
        <div className="grid grid-cols-2 gap-2">
          <FilterInput
            label="Min fitness"
            value={minFitness}
            onChange={(v) => setMinFitness(v === "" ? "" : Number(v))}
            step="0.01"
          />
          <FilterInput
            label="Min ipTM"
            value={minIptm}
            onChange={(v) => setMinIptm(v === "" ? "" : Number(v))}
            step="0.01"
          />
          <FilterInput
            label="Generation"
            value={generation}
            onChange={(v) => setGeneration(v === "" ? "" : Number(v))}
          />
          <label className="flex flex-col gap-1 text-[10px] font-medium text-muted-foreground">
            Lineage
            <select
              value={lineage}
              onChange={(e) => setLineage(e.target.value)}
              className="h-8 rounded-md border bg-background px-2 text-xs"
            >
              <option value="">all</option>
              {lineages.map((l) => (
                <option key={l} value={l}>
                  {l}
                </option>
              ))}
            </select>
          </label>
        </div>
      </div>

      <div className="flex-1 overflow-y-auto p-3">
        {isLoading && (
          <p className="text-center text-xs text-muted-foreground">loading…</p>
        )}
        {error && (
          <p className="text-center text-xs text-red-400">
            {(error as Error).message}
          </p>
        )}
        {!isLoading && filtered.length === 0 && (
          <Card>
            <CardContent className="py-6 text-center text-xs text-muted-foreground">
              No optimized switches match the filters.
            </CardContent>
          </Card>
        )}
        <div className="space-y-3">
          {filtered.map((v) => (
            <VariantCard key={v.variant_id} v={v} />
          ))}
        </div>
      </div>
    </aside>
  );
}

function FilterInput({
  label,
  value,
  onChange,
  step,
}: {
  label: string;
  value: string | number | "";
  onChange: (v: string) => void;
  step?: string;
}) {
  return (
    <label className="flex flex-col gap-1 text-[10px] font-medium text-muted-foreground">
      {label}
      <Input
        type="number"
        step={step}
        value={value as number | ""}
        onChange={(e) => onChange(e.target.value)}
        className="h-8 text-xs"
        placeholder="—"
      />
    </label>
  );
}

function VariantCard({ v }: { v: OptimizedVariant }) {
  const [expanded, setExpanded] = useState(false);
  const [show3D, setShow3D] = useState(false);
  const [state, setState] = useState<StructureState>(() => pickDefaultState(v));
  const openDrawer = useVariantDrawer((s) => s.open);

  const structurePath = structureFor(v, state);

  return (
    <Card>
      <CardHeader className="pb-2">
        <button
          onClick={() => setExpanded((x) => !x)}
          className="flex w-full items-start justify-between gap-2 text-left"
        >
          <div className="flex items-center gap-2">
            {expanded ? (
              <ChevronDown className="mt-0.5 h-3.5 w-3.5 text-muted-foreground" />
            ) : (
              <ChevronRight className="mt-0.5 h-3.5 w-3.5 text-muted-foreground" />
            )}
            <div>
              <CardTitle
                className="cursor-pointer font-mono text-xs hover:underline"
                onClick={(e) => {
                  e.stopPropagation();
                  openDrawer(v.variant_id);
                }}
                title="Open detail drawer"
              >
                {v.variant_id}
              </CardTitle>
              <p className="text-[10px] text-muted-foreground">
                gen {v.generation} · lineage {v.baseline_id ?? "—"}
              </p>
            </div>
          </div>
          <div className="text-right">
            <div className="text-sm font-semibold text-amber-300">
              {fmt(v.fitness, 2)}
            </div>
            <div className="text-[10px] text-muted-foreground">fitness</div>
          </div>
        </button>
      </CardHeader>
      <CardContent className="space-y-2 pt-0">
        <div className="grid grid-cols-3 gap-2 text-[11px]">
          <Metric label="ipTM" value={fmt(v.iptm, 3)} />
          <Metric label="AF2-IG" value={fmt(v.af2_ig, 3)} />
          <Metric label="HEPN Δ" value={fmt(v.hepn_shift_A, 2)} unit="Å" />
          <Metric label="ON dist" value={fmt(v.on_dist_A, 2)} unit="Å" />
          <Metric label="OFF dist" value={fmt(v.off_dist_A, 2)} unit="Å" />
        </div>

        {v.optimized_reasons && v.optimized_reasons.length > 0 && (
          <div className="flex flex-wrap gap-1">
            {v.optimized_reasons.map((r) => (
              <Badge key={r} variant="outline" className="text-[10px]">
                {r}
              </Badge>
            ))}
          </div>
        )}

        {expanded && (
          <>
            <div className="flex items-center justify-between pt-1">
              <div className="flex items-center gap-1 text-[11px]">
                <Box className="h-3 w-3" />
                <select
                  value={state}
                  onChange={(e) => setState(e.target.value as StructureState)}
                  className="h-7 rounded-md border bg-background px-1.5 text-[11px]"
                >
                  <option value="on" disabled={!v.eval_artifacts?.on_structure}>
                    ON
                  </option>
                  <option value="off" disabled={!v.eval_artifacts?.off_structure}>
                    OFF
                  </option>
                  <option
                    value="offtarget"
                    disabled={!v.eval_artifacts?.offtarget_structure}
                  >
                    off-target
                  </option>
                  <option
                    value="optimized"
                    disabled={!v.optimized_artifacts?.structure}
                  >
                    optimized
                  </option>
                </select>
              </div>
              <Button
                size="sm"
                variant={show3D ? "default" : "outline"}
                onClick={() => setShow3D((v) => !v)}
                disabled={!structurePath}
              >
                {show3D ? "Hide 3D" : "Show 3D"}
              </Button>
            </div>
            {show3D && structurePath && (
              <Structure3D
                path={structurePath}
                hepn={v.domain_metadata}
                heightPx={240}
              />
            )}
          </>
        )}
      </CardContent>
    </Card>
  );
}

function Metric({
  label,
  value,
  unit,
}: {
  label: string;
  value: string;
  unit?: string;
}) {
  return (
    <div className="rounded bg-muted/40 px-2 py-1">
      <div className="text-[9px] uppercase text-muted-foreground">{label}</div>
      <div className={cn("font-mono text-[11px]")}>
        {value}
        {unit && <span className="ml-0.5 text-muted-foreground">{unit}</span>}
      </div>
    </div>
  );
}

function fmt(n: number | null | undefined, digits = 2): string {
  if (n == null || Number.isNaN(n)) return "—";
  return Number(n).toFixed(digits);
}

function pickDefaultState(v: OptimizedVariant): StructureState {
  if (v.eval_artifacts?.on_structure) return "on";
  if (v.eval_artifacts?.off_structure) return "off";
  if (v.optimized_artifacts?.structure) return "optimized";
  if (v.eval_artifacts?.offtarget_structure) return "offtarget";
  return "on";
}

function structureFor(
  v: OptimizedVariant,
  state: StructureState,
): string | undefined {
  switch (state) {
    case "on":
      return v.eval_artifacts?.on_structure;
    case "off":
      return v.eval_artifacts?.off_structure;
    case "offtarget":
      return v.eval_artifacts?.offtarget_structure;
    case "optimized":
      return v.optimized_artifacts?.structure;
  }
}
