import { useState, useMemo } from "react";
import { useQuery, useQueryClient } from "@tanstack/react-query";
import { useNavigate } from "react-router-dom";
import { Cpu, Zap, DollarSign, Rocket, RefreshCw } from "lucide-react";

import { api } from "@/lib/api";
import type { Baseline, VastOffer } from "@/lib/types";
import { cn, formatDph } from "@/lib/utils";
import { toast } from "@/components/ui/toast";

import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import {
  Dialog,
  DialogContent,
  DialogHeader,
  DialogTitle,
  DialogDescription,
  DialogFooter,
} from "@/components/ui/dialog";
import { Badge } from "@/components/ui/badge";

interface Props {
  open: boolean;
  onOpenChange: (open: boolean) => void;
  selected: Baseline[];
  onLaunched?: () => void;
}

export function LaunchRunDialog({
  open,
  onOpenChange,
  selected,
  onLaunched,
}: Props) {
  const qc = useQueryClient();
  const navigate = useNavigate();

  const [gpuName, setGpuName] = useState("A100_SXM4");
  const [minVram, setMinVram] = useState(40);
  const [maxDph, setMaxDph] = useState<number | "">(1.5);
  const [offerId, setOfferId] = useState<number | null>(null);

  const [maxGenerations, setMaxGenerations] = useState(12);
  const [variantsPerGen, setVariantsPerGen] = useState(5);
  const [workers, setWorkers] = useState(3);
  const [label, setLabel] = useState("");

  const [launching, setLaunching] = useState(false);

  const offersQ = useQuery({
    enabled: open,
    queryKey: ["vast-offers", { gpuName, minVram, maxDph }],
    queryFn: () =>
      api.listOffers({
        gpu_name: gpuName,
        min_vram_gb: minVram,
        max_dph_total: maxDph === "" ? undefined : maxDph,
        limit: 12,
      }),
  });

  const offers: VastOffer[] = offersQ.data?.offers ?? [];

  const picked = useMemo(
    () => offers.find((o) => o.id === offerId) ?? null,
    [offers, offerId],
  );

  const estimateCostPerHour =
    picked?.dph_total != null ? picked.dph_total * selected.length : null;

  async function handleLaunch() {
    if (!offerId || selected.length === 0) return;
    setLaunching(true);
    try {
      const run = await api.launchRun({
        baseline_ids: selected.map((b) => b.baseline_id),
        crrna_lookup_ids: selected.map((b) => b.crrna_lookup_id),
        offer_id: offerId,
        max_generations: maxGenerations,
        variants_per_gen: variantsPerGen,
        workers,
        label: label.trim() || undefined,
      });
      toast({
        kind: "success",
        title: "Run launched",
        description: `${run.id.slice(0, 8)} • provisioning on Vast.ai`,
      });
      qc.invalidateQueries({ queryKey: ["runs"] });
      onLaunched?.();
      navigate("/runs");
    } catch (e: any) {
      toast({
        kind: "error",
        title: "Launch failed",
        description: String(e?.message ?? e),
      });
    } finally {
      setLaunching(false);
    }
  }

  return (
    <Dialog open={open} onOpenChange={onOpenChange}>
      <DialogContent className="max-h-[85vh] max-w-3xl overflow-y-auto">
        <DialogHeader>
          <DialogTitle className="flex items-center gap-2">
            <Rocket className="h-5 w-5" />
            Launch {selected.length} parallel evolution run
            {selected.length === 1 ? "" : "s"}
          </DialogTitle>
          <DialogDescription>
            Each selected baseline runs on a dedicated Vast.ai A100 VPS using
            the matching offer. Output artifacts will stream back to this
            controller when the run completes.
          </DialogDescription>
        </DialogHeader>

        <section className="space-y-3">
          <div className="flex items-center gap-2 text-sm font-medium">
            <Cpu className="h-4 w-4" /> GPU offer
          </div>
          <div className="grid grid-cols-3 gap-3">
            <LabeledInput
              label="GPU name"
              value={gpuName}
              onChange={(v) => setGpuName(v)}
            />
            <LabeledInput
              label="Min VRAM (GB)"
              value={minVram}
              onChange={(v) => setMinVram(Number(v) || 40)}
              type="number"
            />
            <LabeledInput
              label="Max $/hr"
              value={maxDph === "" ? "" : String(maxDph)}
              onChange={(v) => setMaxDph(v === "" ? "" : Number(v))}
              type="number"
              step="0.01"
            />
          </div>
          <div className="flex items-center justify-between text-xs text-muted-foreground">
            <span>
              {offersQ.isFetching
                ? "searching Vast.ai…"
                : `${offers.length} matching offers`}
            </span>
            <button
              className="inline-flex items-center gap-1 hover:text-foreground"
              onClick={() => offersQ.refetch()}
            >
              <RefreshCw className="h-3 w-3" /> refresh
            </button>
          </div>

          <div className="max-h-64 overflow-y-auto rounded-md border">
            <table className="w-full text-sm">
              <thead className="sticky top-0 bg-muted/80 text-xs uppercase text-muted-foreground">
                <tr>
                  <th className="w-10 p-2" />
                  <th className="p-2 text-left">GPU</th>
                  <th className="p-2 text-left">VRAM</th>
                  <th className="p-2 text-left">CPU / RAM</th>
                  <th className="p-2 text-left">$/h</th>
                  <th className="p-2 text-left">Datacenter</th>
                  <th className="p-2 text-left">Rel.</th>
                </tr>
              </thead>
              <tbody>
                {offers.length === 0 && !offersQ.isFetching && (
                  <tr>
                    <td
                      colSpan={7}
                      className="p-4 text-center text-muted-foreground"
                    >
                      No offers match. Loosen the filters.
                    </td>
                  </tr>
                )}
                {offers.map((o) => (
                  <tr
                    key={o.id}
                    className={cn(
                      "cursor-pointer border-b hover:bg-muted/40",
                      offerId === o.id && "bg-primary/10",
                    )}
                    onClick={() => setOfferId(o.id)}
                  >
                    <td className="p-2">
                      <input
                        type="radio"
                        checked={offerId === o.id}
                        onChange={() => setOfferId(o.id)}
                        className="accent-primary"
                      />
                    </td>
                    <td className="p-2">
                      {o.num_gpus ?? 1}× {o.gpu_name ?? "?"}
                    </td>
                    <td className="p-2">
                      {o.gpu_ram_gb ? `${o.gpu_ram_gb} GB` : "—"}
                    </td>
                    <td className="p-2">
                      {o.cpu_cores ?? "?"}c /{" "}
                      {o.cpu_ram_gb ? `${o.cpu_ram_gb} GB` : "?"}
                    </td>
                    <td className="p-2 font-mono text-xs">
                      {formatDph(o.dph_total)}
                    </td>
                    <td className="p-2 text-xs">{o.datacenter ?? "—"}</td>
                    <td className="p-2 text-xs">
                      {o.reliability != null
                        ? `${(o.reliability * 100).toFixed(1)}%`
                        : "—"}
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </section>

        <section className="space-y-3">
          <div className="flex items-center gap-2 text-sm font-medium">
            <Zap className="h-4 w-4" /> Evolution parameters
          </div>
          <div className="grid grid-cols-4 gap-3">
            <LabeledInput
              label="Max generations"
              type="number"
              value={maxGenerations}
              onChange={(v) => setMaxGenerations(Number(v) || 1)}
            />
            <LabeledInput
              label="Variants / gen"
              type="number"
              value={variantsPerGen}
              onChange={(v) => setVariantsPerGen(Number(v) || 1)}
            />
            <LabeledInput
              label="Workers"
              type="number"
              value={workers}
              onChange={(v) => setWorkers(Number(v) || 1)}
            />
            <LabeledInput
              label="Label (optional)"
              value={label}
              onChange={setLabel}
              placeholder="smoke-test"
            />
          </div>
        </section>

        <section className="flex items-center justify-between rounded-md border bg-muted/30 px-4 py-3 text-sm">
          <div className="flex items-center gap-2">
            <DollarSign className="h-4 w-4" />
            <span>
              {selected.length} runs
              {picked && (
                <>
                  {" @ "}
                  <span className="font-mono">{formatDph(picked.dph_total)}</span>
                </>
              )}
            </span>
          </div>
          {estimateCostPerHour != null && (
            <Badge variant="info">
              ~${estimateCostPerHour.toFixed(2)}/hr fleet cost
            </Badge>
          )}
        </section>

        <DialogFooter>
          <Button variant="ghost" onClick={() => onOpenChange(false)}>
            Cancel
          </Button>
          <Button
            onClick={handleLaunch}
            disabled={!offerId || selected.length === 0 || launching}
          >
            <Rocket className="mr-2 h-4 w-4" />
            {launching ? "Launching…" : "Launch runs"}
          </Button>
        </DialogFooter>
      </DialogContent>
    </Dialog>
  );
}

function LabeledInput({
  label,
  value,
  onChange,
  type = "text",
  step,
  placeholder,
}: {
  label: string;
  value: string | number;
  onChange: (v: string) => void;
  type?: string;
  step?: string;
  placeholder?: string;
}) {
  return (
    <label className="flex flex-col gap-1 text-xs font-medium text-muted-foreground">
      {label}
      <Input
        type={type}
        step={step}
        value={value}
        placeholder={placeholder}
        onChange={(e) => onChange(e.target.value)}
        className="h-9"
      />
    </label>
  );
}
