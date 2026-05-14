import { Component, type ErrorInfo, type ReactNode } from "react";
import { AlertTriangle, RotateCcw } from "lucide-react";

import { Button } from "@/components/ui/button";

interface Props {
  children: ReactNode;
  /** Optional fallback to render; defaults to a centered card. */
  fallback?: (error: Error, reset: () => void) => ReactNode;
}

interface State {
  error: Error | null;
}

/**
 * Phase-3 fix: React used to surface render-time errors as a blank white
 * screen because StrictMode silently swallowed them in dev and crashed the
 * whole tree in prod.  This boundary catches the error, shows a useful
 * message, lets the user click "retry" to remount the subtree, and (in
 * dev) prints the stack to the console so it lands in the Cursor terminal.
 *
 * Mount it at the App root (and re-mount around any subtree that loads
 * untrusted backend data -- the variant drawer, the structure viewer).
 */
export class ErrorBoundary extends Component<Props, State> {
  state: State = { error: null };

  static getDerivedStateFromError(error: Error): State {
    return { error };
  }

  componentDidCatch(error: Error, info: ErrorInfo): void {
    // Keep this visible in the dev console; production builds will still
    // see the user-facing fallback below.
    // eslint-disable-next-line no-console
    console.error("ErrorBoundary caught:", error, info.componentStack);
  }

  reset = () => {
    this.setState({ error: null });
  };

  render(): ReactNode {
    const { error } = this.state;
    if (!error) return this.props.children;

    if (this.props.fallback) return this.props.fallback(error, this.reset);

    return (
      <div className="m-6 max-w-2xl rounded-lg border border-destructive/40 bg-destructive/5 p-6">
        <div className="mb-2 flex items-center gap-2 text-destructive">
          <AlertTriangle className="h-5 w-5" />
          <h2 className="text-sm font-semibold">Something broke in this view</h2>
        </div>
        <p className="mb-3 text-xs text-muted-foreground">
          The CASCADE dashboard encountered a render error. Your data is safe
          and the backend is unaffected; this only impacts the current panel.
        </p>
        <pre className="mb-3 max-h-48 overflow-auto rounded border bg-zinc-950 p-2 text-[11px] text-zinc-300">
          {error.message}
          {error.stack ? `\n\n${error.stack}` : ""}
        </pre>
        <Button variant="outline" size="sm" onClick={this.reset}>
          <RotateCcw className="mr-1 h-3.5 w-3.5" /> Retry
        </Button>
      </div>
    );
  }
}
