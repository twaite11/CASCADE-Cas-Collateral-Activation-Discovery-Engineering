import { clsx, type ClassValue } from "clsx";
import { twMerge } from "tailwind-merge";

export function cn(...inputs: ClassValue[]) {
  return twMerge(clsx(inputs));
}

export function formatDph(dph?: number | null): string {
  if (dph == null) return "—";
  return `$${dph.toFixed(3)}/h`;
}

export function formatRelative(ts?: number | null): string {
  if (!ts) return "—";
  const deltaS = (Date.now() / 1000) - ts;
  if (deltaS < 60) return `${Math.round(deltaS)}s ago`;
  if (deltaS < 3600) return `${Math.round(deltaS / 60)}m ago`;
  if (deltaS < 86400) return `${Math.round(deltaS / 3600)}h ago`;
  return `${Math.round(deltaS / 86400)}d ago`;
}

export function classifyStatus(status: string): "ok" | "warn" | "err" | "info" {
  switch (status) {
    case "completed":
      return "ok";
    case "failed":
      return "err";
    case "cancelled":
      return "warn";
    default:
      return "info";
  }
}
