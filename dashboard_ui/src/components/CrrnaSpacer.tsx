import { cn } from "@/lib/utils";
import { NT_COLORS_TW } from "@/lib/palette";

interface Props {
  repeat?: string | null;
  spacer?: string | null;
}

// F-2: palette imported from `@/lib/palette` so this view and the structure
// viewer never disagree on what "A" or "G" should look like.
const NT_COLORS = NT_COLORS_TW;

/**
 * Monospaced crRNA rendering: 5' — repeat — spacer — 3' with per-nucleotide
 * coloring and a minimap-style track.
 */
export function CrrnaSpacer({ repeat, spacer }: Props) {
  if (!repeat && !spacer) {
    return (
      <p className="text-xs text-muted-foreground">
        No crRNA associated with this baseline.
      </p>
    );
  }
  return (
    <div className="space-y-2">
      <div className="flex flex-wrap items-baseline gap-2 font-mono text-[11px]">
        <span className="text-muted-foreground">5′</span>
        {repeat && (
          <span className="rounded bg-muted/50 px-1 py-0.5" title="Direct repeat">
            {colorize(repeat)}
          </span>
        )}
        {spacer && (
          <span
            className="rounded bg-primary/10 px-1 py-0.5"
            title="Spacer (target match)"
          >
            {colorize(spacer)}
          </span>
        )}
        <span className="text-muted-foreground">3′</span>
      </div>
      <div className="flex gap-0.5 text-[9px]">
        {repeat && (
          <TrackLabel count={repeat.length} tone="muted" label="DR" />
        )}
        {spacer && (
          <TrackLabel count={spacer.length} tone="primary" label="spacer" />
        )}
      </div>
    </div>
  );
}

function colorize(seq: string) {
  return seq.split("").map((nt, i) => (
    <span key={i} className={cn(NT_COLORS[nt.toUpperCase()] ?? "")}>
      {nt}
    </span>
  ));
}

function TrackLabel({
  count,
  tone,
  label,
}: {
  count: number;
  tone: "muted" | "primary";
  label: string;
}) {
  return (
    <div
      className={cn(
        "flex-1 rounded border-b-2 px-1 text-[9px]",
        tone === "primary"
          ? "border-primary/60 text-primary"
          : "border-muted-foreground/40 text-muted-foreground",
      )}
      style={{ flexGrow: count }}
    >
      {label} ({count} nt)
    </div>
  );
}
