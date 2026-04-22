import { useMemo, useState } from "react";
import { useQuery } from "@tanstack/react-query";
import { Search, FlaskConical, ShieldCheck, Box, Rocket } from "lucide-react";

import { api } from "@/lib/api";
import type { Baseline } from "@/lib/types";
import { cn } from "@/lib/utils";

import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { Switch } from "@/components/ui/switch";
import { Badge } from "@/components/ui/badge";
import {
  Card,
  CardContent,
  CardHeader,
  CardTitle,
} from "@/components/ui/card";
import { LaunchRunDialog } from "@/components/LaunchRunDialog";

export function BaselinesTable() {
  const [search, setSearch] = useState("");
  const [validatedOnly, setValidatedOnly] = useState(false);
  const [withStructureOnly, setWithStructureOnly] = useState(false);
  const [selected, setSelected] = useState<Set<string>>(new Set());
  const [launchOpen, setLaunchOpen] = useState(false);

  const { data, isLoading, error } = useQuery({
    queryKey: [
      "baselines",
      { search, validatedOnly, withStructureOnly },
    ],
    queryFn: () =>
      api.listBaselines({
        search: search.trim() || undefined,
        validated_only: validatedOnly || undefined,
        with_structure_only: withStructureOnly || undefined,
        limit: 500,
      }),
  });

  const rows = data?.rows ?? [];
  const selectedRows = useMemo(
    () => rows.filter((r) => selected.has(r.baseline_id)),
    [rows, selected],
  );

  function toggle(id: string) {
    setSelected((prev) => {
      const next = new Set(prev);
      if (next.has(id)) next.delete(id);
      else next.add(id);
      return next;
    });
  }

  function toggleAllVisible() {
    setSelected((prev) => {
      const next = new Set(prev);
      const allChecked = rows.every((r) => next.has(r.baseline_id));
      if (allChecked) rows.forEach((r) => next.delete(r.baseline_id));
      else rows.forEach((r) => next.add(r.baseline_id));
      return next;
    });
  }

  return (
    <div className="space-y-4">
      <Card>
        <CardHeader className="flex-row items-center justify-between space-y-0 pb-4">
          <div>
            <CardTitle className="flex items-center gap-2 text-base">
              <FlaskConical className="h-4 w-4" />
              Enzyme baselines
              {data && (
                <span className="text-xs font-normal text-muted-foreground">
                  {data.total} total
                </span>
              )}
            </CardTitle>
          </div>
          <Button
            disabled={selected.size === 0}
            onClick={() => setLaunchOpen(true)}
          >
            <Rocket className="mr-2 h-4 w-4" />
            Launch {selected.size || ""} parallel run
            {selected.size === 1 ? "" : "s"}
          </Button>
        </CardHeader>
        <CardContent className="space-y-3">
          <div className="flex flex-wrap items-center gap-4">
            <div className="relative min-w-[260px] flex-1">
              <Search className="absolute left-3 top-1/2 h-4 w-4 -translate-y-1/2 text-muted-foreground" />
              <Input
                placeholder="Search enzyme IDs, SRA accession…"
                className="pl-9"
                value={search}
                onChange={(e) => setSearch(e.target.value)}
              />
            </div>
            <label className="flex items-center gap-2 text-sm">
              <Switch
                checked={validatedOnly}
                onCheckedChange={setValidatedOnly}
              />
              <ShieldCheck className="h-4 w-4" /> Validated only
            </label>
            <label className="flex items-center gap-2 text-sm">
              <Switch
                checked={withStructureOnly}
                onCheckedChange={setWithStructureOnly}
              />
              <Box className="h-4 w-4" /> Has Phase 1 structure
            </label>
          </div>

          <div className="overflow-hidden rounded-md border">
            <table className="w-full text-sm">
              <thead className="border-b bg-muted/50 text-xs uppercase text-muted-foreground">
                <tr>
                  <th className="w-10 p-2 text-left">
                    <input
                      type="checkbox"
                      checked={
                        rows.length > 0 &&
                        rows.every((r) => selected.has(r.baseline_id))
                      }
                      onChange={toggleAllVisible}
                      className="h-4 w-4 accent-primary"
                      aria-label="Select all visible"
                    />
                  </th>
                  <th className="p-2 text-left">Baseline ID</th>
                  <th className="p-2 text-left">Subtype</th>
                  <th className="p-2 text-left">Length</th>
                  <th className="p-2 text-left">HEPN1</th>
                  <th className="p-2 text-left">HEPN2</th>
                  <th className="p-2 text-left">crRNA</th>
                  <th className="p-2 text-left">Flags</th>
                </tr>
              </thead>
              <tbody>
                {isLoading && (
                  <tr>
                    <td colSpan={8} className="p-6 text-center text-muted-foreground">
                      Loading baselines…
                    </td>
                  </tr>
                )}
                {error && (
                  <tr>
                    <td colSpan={8} className="p-6 text-center text-red-400">
                      Failed to load baselines: {(error as Error).message}
                    </td>
                  </tr>
                )}
                {!isLoading && rows.length === 0 && (
                  <tr>
                    <td colSpan={8} className="p-6 text-center text-muted-foreground">
                      No baselines match your filters.
                    </td>
                  </tr>
                )}
                {rows.map((b) => (
                  <BaselineRow
                    key={b.baseline_id}
                    row={b}
                    checked={selected.has(b.baseline_id)}
                    onToggle={() => toggle(b.baseline_id)}
                  />
                ))}
              </tbody>
            </table>
          </div>
        </CardContent>
      </Card>

      <LaunchRunDialog
        open={launchOpen}
        onOpenChange={setLaunchOpen}
        selected={selectedRows}
        onLaunched={() => {
          setSelected(new Set());
          setLaunchOpen(false);
        }}
      />
    </div>
  );
}

function BaselineRow({
  row,
  checked,
  onToggle,
}: {
  row: Baseline;
  checked: boolean;
  onToggle: () => void;
}) {
  return (
    <tr
      className={cn(
        "border-b transition-colors hover:bg-muted/40",
        checked && "bg-primary/5",
      )}
    >
      <td className="p-2">
        <input
          type="checkbox"
          checked={checked}
          onChange={onToggle}
          className="h-4 w-4 accent-primary"
          aria-label={`Select ${row.baseline_id}`}
        />
      </td>
      <td className="p-2 font-mono text-xs">
        <button
          onClick={onToggle}
          className="text-left hover:underline"
          title={row.sra_accession ?? undefined}
        >
          {row.baseline_id}
        </button>
      </td>
      <td className="p-2">{row.subtype ?? "—"}</td>
      <td className="p-2">{row.sequence_length ?? "—"}</td>
      <td className="p-2">
        {row.hepn1_start != null && row.hepn1_end != null
          ? `${row.hepn1_start}–${row.hepn1_end}`
          : "—"}
      </td>
      <td className="p-2">
        {row.hepn2_start != null && row.hepn2_end != null
          ? `${row.hepn2_start}–${row.hepn2_end}`
          : "—"}
      </td>
      <td className="p-2">
        {row.crrna_spacer ? (
          <code
            className="truncate rounded bg-muted px-1 text-[11px]"
            title={row.crrna_spacer}
          >
            {row.crrna_spacer.slice(0, 14)}
            {row.crrna_spacer.length > 14 ? "…" : ""}
          </code>
        ) : (
          "—"
        )}
      </td>
      <td className="p-2">
        <div className="flex flex-wrap gap-1">
          {row.validated && <Badge variant="ok">validated</Badge>}
          {row.has_phase1_structure && <Badge variant="info">phase1</Badge>}
          {row.status === "elite" && <Badge variant="warn">elite</Badge>}
        </div>
      </td>
    </tr>
  );
}
