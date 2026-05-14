import { useState, useMemo, useEffect } from "react";
import { useQuery, useQueryClient } from "@tanstack/react-query";
import { useNavigate } from "react-router-dom";
import { Cpu, Zap, DollarSign, Rocket, RefreshCw, Split } from "lucide-react";

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
  /** Pooled launch or baseline #0 in fan-out with "same offer" lock. */
  const [offerId, setOfferId] = useState<number | null>(null);
  /** Fan-out only: baseline_id -> vast offer ask id (separate GPUs / racks). */
  const [offerByBaselineId, setOfferByBaselineId] = useState<Record<string, number>>({});
  /** When true (default ON for fan-out), each parallel POST uses its own selected offer ID. */
  const [distinctOffersPerBaseline, setDistinctOffersPerBaseline] = useState(true);

  const [maxGenerations, setMaxGenerations] = useState(12);
  const [variantsPerGen, setVariantsPerGen] = useState(5);
  const [workers, setWorkers] = useState(3);
  const [label, setLabel] = useState("");

  // Option B fan-out: when ON (and >1 baseline selected), each baseline
  // gets its OWN /api/runs POST -> its own Vast.ai VPS -> its own log
  // stream.  When OFF, the legacy single-VPS multi-worker behaviour is
  // used (workers=N, one VPS).  Defaults ON for >=2 baselines so the
  // dialog's "dedicated VPS per baseline" copy is actually true.
  const [fanOut, setFanOut] = useState(true);
  const [launching, setLaunching] = useState(false);
  const [launchProgress, setLaunchProgress] = useState<{ done: number; total: number } | null>(null);

  const offersQ = useQuery({
    enabled: open,
    queryKey: ["vast-offers", { gpuName, minVram, maxDph, wideList: fanOut }],
    queryFn: () =>
      api.listOffers({
        gpu_name: gpuName,
        min_vram_gb: minVram,
        max_dph_total: maxDph === "" ? undefined : maxDph,
        limit: fanOut ? 36 : 12,
      }),
  });

  const offers: VastOffer[] = offersQ.data?.offers ?? [];

  const picked = useMemo(
    () => offers.find((o) => o.id === offerId) ?? null,
    [offers, offerId],
  );

  // In fan-out mode each baseline is its own Vast.ai VPS, so the fleet
  // cost is N x dph.  In single-VPS mode there's just one VPS.
  const willFanOut = fanOut && selected.length > 1;
  const fleetSize = willFanOut ? selected.length : 1;

  const selectedKey = useMemo(
    () => selected.map((s) => s.baseline_id).join("|"),
    [selected],
  );

  const offerIdFingerprint = useMemo(
    () => offers.map((o) => o.id).join(","),
    [offers],
  );

  // Fan-out + distinct offers: stagger defaults across rows so parallel
  // rentals target different GPU slots instead of cloning the same ask id three times.
  useEffect(() => {
    if (!open || offers.length === 0) return;
    if (!willFanOut || !distinctOffersPerBaseline || selected.length === 0)
      return;
    setOfferByBaselineId((prev) => {
      const next: Record<string, number> = {};
      for (let idx = 0; idx < selected.length; idx++) {
        const b = selected[idx];
        const old = prev[b.baseline_id];
        const retains =
          typeof old === "number" && offers.some((o) => o.id === old);
        if (retains) {
          next[b.baseline_id] = old!;
        } else {
          next[b.baseline_id] =
            offers[idx % offers.length]?.id ?? offers[0]!.id;
        }
      }
      return next;
    });
  }, [
    open,
    willFanOut,
    distinctOffersPerBaseline,
    selectedKey,
    offerIdFingerprint,
  ]);

  useEffect(() => {
    if (!open || offers.length === 0) return;
    if (willFanOut && distinctOffersPerBaseline) return;
    const currentValid =
      offerId != null && offers.some((o) => o.id === offerId);
    if (!currentValid) {
      setOfferId(offers[0]?.id ?? null);
    }
  }, [open, willFanOut, distinctOffersPerBaseline, offerIdFingerprint, offerId]);

  const estimateCostPerHour = useMemo(() => {
    if (offers.length === 0 || selected.length === 0) return null;
    if (willFanOut && distinctOffersPerBaseline && selected.length > 1) {
      let sum = 0;
      for (const b of selected) {
        const oid = offerByBaselineId[b.baseline_id];
        const o =
          oid != null ? offers.find((x) => x.id === oid) : undefined;
        if (o?.dph_total == null) return null;
        sum += o.dph_total;
      }
      return sum;
    }
    if (picked?.dph_total == null) return null;
    return picked.dph_total * fleetSize;
  }, [
    offers,
    selected,
    willFanOut,
    distinctOffersPerBaseline,
    offerByBaselineId,
    picked,
    fleetSize,
  ]);

  const launchOffersReady =
    selected.length > 0 &&
    (willFanOut && distinctOffersPerBaseline
      ? selected.every((b) => typeof offerByBaselineId[b.baseline_id] === "number")
      : offerId != null);

  async function handleLaunch() {
    if (!launchOffersReady) return;
    setLaunching(true);
    setLaunchProgress(null);
    try {
      if (willFanOut) {
        // Option B: one POST per baseline, in parallel, each with
        // workers=1 since each VPS will only see one lineage.  The
        // controller assigns a fresh run_id per call so the runs.db
        // rows, log channels, and artifacts dirs are all isolated.
        setLaunchProgress({ done: 0, total: selected.length });
        let done = 0;
        const results = await Promise.allSettled(
          selected.map((b) =>
            api
              .launchRun({
                baseline_ids: [b.baseline_id],
                crrna_lookup_ids: [b.crrna_lookup_id],
                offer_id:
                  distinctOffersPerBaseline
                    ? offerByBaselineId[b.baseline_id]!
                    : offerId!,
                max_generations: maxGenerations,
                variants_per_gen: variantsPerGen,
                workers: 1,
                label: label.trim()
                  ? `${label.trim()} • ${b.baseline_id}`
                  : b.baseline_id,
              })
              .then((r) => {
                done += 1;
                setLaunchProgress({ done, total: selected.length });
                return r;
              }),
          ),
        );
        const ok = results.filter((r) => r.status === "fulfilled").length;
        const fail = results.length - ok;
        if (ok > 0) {
          toast({
            kind: fail === 0 ? "success" : "info",
            title: `Launched ${ok}/${selected.length} runs`,
            description:
              fail === 0
                ? "All baselines provisioning in parallel on Vast.ai"
                : `${fail} launch(es) failed — see /runs for details`,
          });
        }
        if (fail > 0 && ok === 0) {
          // Surface the first error so the user gets actionable info.
          const firstFail = results.find((r) => r.status === "rejected");
          const reason = firstFail?.status === "rejected" ? String(firstFail.reason?.message ?? firstFail.reason) : "unknown";
          throw new Error(reason);
        }
      } else {
        // Legacy single-VPS mode (Option A): N baselines share one VPS
        // and the orchestrator parallelises them via N worker processes.
        const run = await api.launchRun({
          baseline_ids: selected.map((b) => b.baseline_id),
          crrna_lookup_ids: selected.map((b) => b.crrna_lookup_id),
          offer_id: offerId!,
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
      }
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
      setLaunchProgress(null);
    }
  }

  return (
    <Dialog open={open} onOpenChange={onOpenChange}>
      <DialogContent className="max-h-[85vh] max-w-3xl overflow-y-auto">
        <DialogHeader>
          <DialogTitle className="flex items-center gap-2">
            <Rocket className="h-5 w-5" />
            Launch {willFanOut ? selected.length : 1} evolution run
            {willFanOut && selected.length !== 1 ? "s" : ""}
            {willFanOut ? " (fan-out)" : ""}
          </DialogTitle>
          <DialogDescription>
            {willFanOut
              ? `Each of the ${selected.length} selected baselines will provision its own Vast.ai VPS and run with workers=1. Three independent log streams, three independent artifact dirs, true wall-clock parallelism.`
              : selected.length > 1
                ? `All ${selected.length} baselines will share a single Vast.ai VPS and run as ${workers} workers on that host (orchestrator multiprocessing).`
                : "Single baseline will run on the selected Vast.ai VPS."}
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
                {offers.map((o) => {
                  const useBrowseRadios =
                    !willFanOut || !distinctOffersPerBaseline;
                  return (
                  <tr
                    key={o.id}
                    className={cn(
                      useBrowseRadios && "cursor-pointer hover:bg-muted/40",
                      "border-b",
                      offerId === o.id && useBrowseRadios && "bg-primary/10",
                    )}
                    onClick={
                      useBrowseRadios ? () => setOfferId(o.id) : undefined
                    }
                  >
                    <td className="p-2">
                      {useBrowseRadios ? (
                      <input
                        type="radio"
                        checked={offerId === o.id}
                        onChange={() => setOfferId(o.id)}
                        className="accent-primary"
                      />
                      ) : (
                      <span className="inline-block h-4 w-4 px-px text-muted-foreground"
                        aria-hidden />
                      )}
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
                );
                })}
              </tbody>
            </table>
            {willFanOut && distinctOffersPerBaseline && (
              <p className="text-xs text-muted-foreground px-2 pt-1">
                With per-enzyme offers enabled, the GPU table above is for
                comparison only — assign each VPS using the selectors below.
              </p>
            )}
          </div>
        </section>

        {selected.length > 1 && (
          <section className="space-y-2 rounded-md border border-primary/30 bg-primary/5 px-4 py-3">
            <label className="flex cursor-pointer items-start gap-3 text-sm">
              <input
                type="checkbox"
                checked={fanOut}
                onChange={(e) => setFanOut(e.target.checked)}
                className="mt-0.5 h-4 w-4 accent-primary"
              />
              <div className="space-y-1">
                <div className="flex items-center gap-2 font-medium">
                  <Split className="h-4 w-4" />
                  Fan out: one Vast.ai VPS per baseline ({selected.length}×)
                </div>
                <p className="text-xs text-muted-foreground">
                  ON (recommended): {selected.length} independent VPS instances,
                  each running one lineage, full wall-clock parallelism.
                  Costs {selected.length}× per hour but finishes in ~1× wall time.
                  <br />
                  OFF: legacy mode — one VPS shared across {selected.length} workers,
                  cheaper but slower (GPU lock serialises Protenix calls).
                </p>
              </div>
            </label>

            {willFanOut && (
              <>
                <label className="ml-10 flex cursor-pointer items-start gap-3 text-xs">
                  <input
                    type="checkbox"
                    checked={distinctOffersPerBaseline}
                    onChange={(e) => {
                      const on = e.target.checked;
                      setDistinctOffersPerBaseline(on);
                      if (!on && offers.length && selected.length) {
                        const seed =
                          offerByBaselineId[selected[0]!.baseline_id] ??
                          offerId ??
                          offers[0]?.id ??
                          null;
                        if (seed != null) setOfferId(seed);
                      }
                    }}
                    disabled={launching}
                    className="mt-1 h-3.5 w-3.5 accent-primary"
                  />
                  <div>
                    <span className="font-medium">
                      Dedicated Vast.ai offer ID per enzyme
                    </span>
                    <p className="text-muted-foreground">
                      Default ON: each POST sends a different Vast.ai offer ID
                      (chosen below)—less contention than cloning one ask multiple
                      times. Turn off only if every run must share one listing.
                    </p>
                  </div>
                </label>

                {distinctOffersPerBaseline && offers.length > 0 && (
                  <div className="ml-10 mt-4 space-y-2 rounded-md border border-border bg-muted/20 p-3">
                    <div className="text-xs font-medium uppercase text-muted-foreground">
                      Assign offer per enzyme
                    </div>
                    <div className="max-h-48 space-y-2 overflow-y-auto">
                      {selected.map((b) => (
                        <div
                          key={b.baseline_id}
                          className="flex flex-col gap-1 sm:flex-row sm:items-center sm:gap-3"
                        >
                          <label
                            className="shrink-0 font-mono text-[11px] text-muted-foreground sm:w-[30%]"
                            title={b.baseline_id}
                          >
                            {baselineShort(b.baseline_id)}
                          </label>
                          <select
                            value={offerByBaselineId[b.baseline_id] ?? ""}
                            onChange={(ev) =>
                              setOfferByBaselineId((prev) => ({
                                ...prev,
                                [b.baseline_id]: Number(ev.target.value),
                              }))
                            }
                            className="h-9 w-full rounded-md border border-input bg-background px-2 text-sm"
                          >
                            {offers.map((o) => (
                              <option key={o.id} value={o.id}>
                                #{o.id} · {formatDph(o.dph_total)} ·{" "}
                                {o.datacenter ?? "?"} ({o.gpu_name ?? "GPU"})
                              </option>
                            ))}
                          </select>
                        </div>
                      ))}
                    </div>
                  </div>
                )}
              </>
            )}
          </section>
        )}

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
              label={willFanOut ? "Workers (per VPS)" : "Workers"}
              type="number"
              value={willFanOut ? 1 : workers}
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
              {fleetSize} VPS instance{fleetSize === 1 ? "" : "s"}
              {willFanOut &&
              distinctOffersPerBaseline &&
              selected.length > 1 ? (
                <span className="text-muted-foreground">
                  {" "}(each offer chosen below)
                </span>
              ) : (
                picked && (
                  <>
                    {" @ "}
                    <span className="font-mono">{formatDph(picked.dph_total)}</span>
                    {fleetSize > 1 && <span> each</span>}
                  </>
                )
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
          <Button variant="ghost" onClick={() => onOpenChange(false)} disabled={launching}>
            Cancel
          </Button>
          <Button
            onClick={handleLaunch}
            disabled={!launchOffersReady || launching}
          >
            <Rocket className="mr-2 h-4 w-4" />
            {launching
              ? launchProgress
                ? `Launching ${launchProgress.done}/${launchProgress.total}…`
                : "Launching…"
              : willFanOut
                ? `Launch ${selected.length} runs`
                : "Launch run"}
          </Button>
        </DialogFooter>
      </DialogContent>
    </Dialog>
  );
}

function baselineShort(full: string) {
  if (full.length <= 48) return full;
  return `${full.slice(0, 22)}…${full.slice(-20)}`;
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
