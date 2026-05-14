import { useEffect, useRef } from "react";
import { Terminal } from "@xterm/xterm";
import { FitAddon } from "@xterm/addon-fit";
import "@xterm/xterm/css/xterm.css";

import { openLogStream } from "@/lib/api";
import { PALETTE_HEX } from "@/lib/palette";

interface Props {
  runId: string;
  initialLines?: string[];
  heightPx?: number;
}

/**
 * xterm.js terminal bound to `/api/runs/{id}/logs` WebSocket stream.
 *
 * C-16 fix: the previous implementation seeded the buffer with
 * `initialLines` *exactly once* on mount and then `eslint-disable`d the
 * dependency array.  When the parent component swapped between two runs,
 * the prop changed but the seed never re-ran -- the terminal kept the old
 * run's history with new live lines interleaved.
 *
 * We now key the entire useEffect on `runId` AND re-seed using a ref so
 * the most-recent `initialLines` is always reflected without triggering
 * the full xterm remount.  We also key on the buffer's content hash so
 * an HTTP snapshot refresh after the user re-opens the tab actually
 * back-fills lines the WS missed during disconnect.
 */
export function LogTerminal({ runId, initialLines = [], heightPx = 320 }: Props) {
  const hostRef = useRef<HTMLDivElement | null>(null);
  const termRef = useRef<Terminal | null>(null);
  const fitRef = useRef<FitAddon | null>(null);
  // Track the seed we've already written so we don't duplicate lines when
  // the parent refetches the HTTP snapshot.
  const seedSignatureRef = useRef<string>("");

  // Re-seed whenever the run changes or the snapshot grows past what we've
  // already written.  Cheap content fingerprint = "len|first|last".
  useEffect(() => {
    if (!termRef.current || initialLines.length === 0) return;
    const sig = `${initialLines.length}|${initialLines[0]}|${
      initialLines[initialLines.length - 1]
    }`;
    if (sig === seedSignatureRef.current) return;
    // We've grown; only write the *new tail*, not the whole buffer again.
    const oldLen = parseInt(seedSignatureRef.current.split("|")[0] ?? "0", 10);
    const tail = initialLines.slice(Number.isFinite(oldLen) ? oldLen : 0);
    for (const line of tail) termRef.current.writeln(line);
    seedSignatureRef.current = sig;
  }, [initialLines]);

  useEffect(() => {
    if (!hostRef.current) return;

    const term = new Terminal({
      convertEol: true,
      fontFamily:
        "JetBrains Mono, Menlo, Consolas, 'Courier New', monospace",
      fontSize: 12,
      theme: {
        background: PALETTE_HEX.structureBg,
        foreground: "#e4e4e7",
        cursor: PALETTE_HEX.rna,
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
    // Reset seed tracking and write the current snapshot for this run.
    seedSignatureRef.current = "";
    if (initialLines.length > 0) {
      for (const line of initialLines) term.writeln(line);
      seedSignatureRef.current = `${initialLines.length}|${initialLines[0]}|${
        initialLines[initialLines.length - 1]
      }`;
    }

    const close = openLogStream(
      runId,
      (line) => term.writeln(line),
      (reason) =>
        term.writeln(
          `\x1b[2m[connection closed${reason ? `: ${reason}` : ""}]\x1b[0m`,
        ),
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
    // Only re-mount when runId changes; initialLines re-seed via the
    // dedicated effect above.
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
