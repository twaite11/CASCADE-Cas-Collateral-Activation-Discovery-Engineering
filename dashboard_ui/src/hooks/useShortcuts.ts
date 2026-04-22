import { useEffect, useRef } from "react";
import { useNavigate } from "react-router-dom";

import { useHelp } from "@/components/HelpPanel";
import { useVariantDrawer } from "@/stores/variantDrawer";

const NAV_MAP: Record<string, string> = {
  r: "/runs",
  b: "/baselines",
  v: "/variants",
  o: "/optimized",
  p: "/production",
};

/**
 * Global key bindings for the dashboard. Press `?` to open the help panel,
 * `g r/b/v/o/p` to jump between tabs, `/` to focus the nearest search input,
 * and `Esc` to close the help panel or the variant drawer.
 */
export function useShortcuts() {
  const navigate = useNavigate();
  const toggleHelp = useHelp((s) => s.toggle);
  const setHelpOpen = useHelp((s) => s.set);
  const closeDrawer = useVariantDrawer((s) => s.close);
  const pendingG = useRef<number | null>(null);

  useEffect(() => {
    function onKey(e: KeyboardEvent) {
      // Never hijack keystrokes while the user is typing in a form field.
      const target = e.target as HTMLElement | null;
      if (
        target &&
        (target.tagName === "INPUT" ||
          target.tagName === "TEXTAREA" ||
          target.tagName === "SELECT" ||
          target.isContentEditable)
      ) {
        if (e.key === "Escape") target.blur();
        return;
      }

      if (e.key === "?" || (e.shiftKey && e.key === "/")) {
        e.preventDefault();
        toggleHelp();
        return;
      }
      if (e.key === "Escape") {
        setHelpOpen(false);
        closeDrawer();
        return;
      }
      if (e.key === "/") {
        const search = document.querySelector<HTMLInputElement>(
          "input[placeholder*='Search' i]",
        );
        if (search) {
          e.preventDefault();
          search.focus();
        }
        return;
      }

      if (e.key === "g") {
        pendingG.current = Date.now();
        setTimeout(() => {
          pendingG.current = null;
        }, 800);
        return;
      }

      if (pendingG.current && Date.now() - pendingG.current < 800) {
        const route = NAV_MAP[e.key];
        if (route) {
          e.preventDefault();
          navigate(route);
          pendingG.current = null;
        }
      }
    }

    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [navigate, toggleHelp, setHelpOpen, closeDrawer]);
}
