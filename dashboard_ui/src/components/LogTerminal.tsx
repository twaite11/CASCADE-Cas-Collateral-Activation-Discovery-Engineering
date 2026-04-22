import { useEffect, useRef } from "react";
import { Terminal } from "@xterm/xterm";
import { FitAddon } from "@xterm/addon-fit";
import "@xterm/xterm/css/xterm.css";

import { openLogStream } from "@/lib/api";

interface Props {
  runId: string;
  initialLines?: string[];
  heightPx?: number;
}

/**
 * xterm.js terminal bound to `/api/runs/{id}/logs` WebSocket stream. Seeds
 * with `initialLines` (from HTTP snapshot), then streams live lines.
 */
export function LogTerminal({ runId, initialLines = [], heightPx = 320 }: Props) {
  const hostRef = useRef<HTMLDivElement | null>(null);
  const termRef = useRef<Terminal | null>(null);
  const fitRef = useRef<FitAddon | null>(null);

  useEffect(() => {
    if (!hostRef.current) return;

    const term = new Terminal({
      convertEol: true,
      fontFamily:
        "JetBrains Mono, Menlo, Consolas, 'Courier New', monospace",
      fontSize: 12,
      theme: {
        background: "#09090b",
        foreground: "#e4e4e7",
        cursor: "#38bdf8",
        selectionBackground: "#1e40af",
      },
      scrollback: 10_000,
      disableStdin: true,
    });
    const fit = new FitAddon();
    term.loadAddon(fit);
    term.open(hostRef.current);
    try {
      fit.fit();
    } catch {
      /* noop */
    }

    termRef.current = term;
    fitRef.current = fit;

    for (const line of initialLines) {
      term.writeln(line);
    }

    const close = openLogStream(
      runId,
      (line) => term.writeln(line),
      () => term.writeln("\x1b[2m[connection closed]\x1b[0m"),
    );

    const ro = new ResizeObserver(() => {
      try {
        fit.fit();
      } catch {
        /* noop */
      }
    });
    ro.observe(hostRef.current);

    return () => {
      ro.disconnect();
      close();
      term.dispose();
      termRef.current = null;
      fitRef.current = null;
    };
    // Only re-mount when runId changes — initialLines re-renders are fine to ignore.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [runId]);

  return (
    <div
      ref={hostRef}
      className="overflow-hidden rounded-md border bg-zinc-950"
      style={{ height: heightPx }}
    />
  );
}
