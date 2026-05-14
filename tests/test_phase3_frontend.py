"""Phase 3 frontend regression tests.

The TS surface is verified by `npm run build` (tsc + vite) in CI; this
suite only checks that the source files we just added are present, the
shared palette is referenced by every component that needs HEPN colors,
and the stub routes are no longer pure stubs.

Lightweight by design: full DOM tests would require a Node test runner
which we don't yet have wired up.
"""
from __future__ import annotations

from pathlib import Path

import pytest

UI = Path(__file__).resolve().parent.parent / "dashboard_ui" / "src"


def _read(rel: str) -> str:
    return (UI / rel).read_text(encoding="utf-8")


# ---------------------------------------------------------------------------
# C-10 / C-11 normalizer
# ---------------------------------------------------------------------------
class TestC10_C11_Adapter:
    def test_adapter_module_exists(self):
        assert (UI / "lib" / "adapters.ts").exists()

    def test_normalize_variant_exported(self):
        src = _read("lib/adapters.ts")
        assert "export function normalizeVariant" in src
        assert "NormalizedVariant" in src
        assert "NormalizedHepn" in src
        assert "NormalizedCrrna" in src

    def test_adapter_consolidates_known_sources(self):
        src = _read("lib/adapters.ts")
        # crRNA alias resolver checks for both snake_case and camelCase keys.
        for alias in ("crrna_repeat", "crRNA_repeat", "crrna_spacer", "crRNA_spacer"):
            assert alias in src
        # HEPN keys
        for key in ("hepn1_start", "hepn1_end", "hepn2_start", "hepn2_end"):
            assert key in src

    def test_drawer_uses_normalizer(self):
        src = _read("components/VariantDetailDrawer.tsx")
        assert "normalizeVariant" in src
        assert "NormalizedVariant" in src
        # Strip line comments so we don't false-match on the explanatory
        # block describing the old buggy pattern.
        code = "\n".join(
            line for line in src.splitlines() if not line.lstrip().startswith("//")
        )
        # The hand-rolled `(v.domain_metadata as any)?.crrna_repeat` pattern
        # must be gone -- it's the exact bug this phase fixes.
        assert "(v.domain_metadata as any)?.crrna_repeat" not in code
        assert "(v.catalog_metadata as any)?.crrna_repeat" not in code


# ---------------------------------------------------------------------------
# F-2 palette parity
# ---------------------------------------------------------------------------
class TestF2_PaletteParity:
    def test_palette_module_exists(self):
        assert (UI / "lib" / "palette.ts").exists()
        src = _read("lib/palette.ts")
        assert "PALETTE_HEX" in src
        assert "HEPN_LEGEND_ITEMS" in src

    def test_structure3d_uses_shared_palette(self):
        src = _read("components/Structure3D.tsx")
        assert 'from "@/lib/palette"' in src
        assert "PALETTE_HEX" in src
        # Inline hex codes for HEPN should no longer appear.
        assert '"#f97316"' not in src
        assert '"#ef4444"' not in src

    def test_crrna_spacer_uses_shared_palette(self):
        src = _read("components/CrrnaSpacer.tsx")
        assert 'from "@/lib/palette"' in src
        assert "NT_COLORS_TW" in src


# ---------------------------------------------------------------------------
# C-16 LogTerminal seed buffer
# ---------------------------------------------------------------------------
class TestC16_LogTerminalSeed:
    def test_seeds_via_signature_so_resnaps_dont_duplicate(self):
        src = _read("components/LogTerminal.tsx")
        # The fix introduces a content fingerprint so re-snapshotting only
        # writes the tail; verify the markers.
        assert "seedSignatureRef" in src
        assert "writeln" in src
        # The old "eslint-disable" + single-seed-on-mount comment should be
        # replaced with the dedicated effect describing the fix.
        assert "C-16 fix" in src


# ---------------------------------------------------------------------------
# C-2 client auth
# ---------------------------------------------------------------------------
class TestC2_ClientAuth:
    def test_api_attaches_key_header(self):
        src = _read("lib/api.ts")
        assert "X-Cascade-Api-Key" in src
        assert "getApiKey" in src
        assert "withAuth" in src

    def test_ws_attaches_key_query_param(self):
        src = _read("lib/api.ts")
        # Browsers can't set custom WS headers; must use ?api_key= query.
        assert "api_key=" in src


# ---------------------------------------------------------------------------
# F-1 stub routes filled in
# ---------------------------------------------------------------------------
class TestF1_StubRoutesFilled:
    @pytest.mark.parametrize(
        "file_path,banned",
        [
            ("routes/VariantsPage.tsx", "Coming soon"),
            ("routes/OptimizedPage.tsx", "wires up in a subsequent commit"),
            ("routes/ProductionPage.tsx", "Coming soon"),
        ],
    )
    def test_no_placeholder_text(self, file_path, banned):
        src = _read(file_path)
        assert banned not in src, f"{file_path} still has placeholder text {banned!r}"

    def test_routes_actually_fetch_data(self):
        for f in (
            "routes/VariantsPage.tsx",
            "routes/OptimizedPage.tsx",
            "routes/ProductionPage.tsx",
        ):
            src = _read(f)
            assert "useQuery" in src, f"{f} should fetch data via React Query"


# ---------------------------------------------------------------------------
# Error boundary wired in
# ---------------------------------------------------------------------------
class TestErrorBoundary:
    def test_boundary_component_exists(self):
        assert (UI / "components" / "ErrorBoundary.tsx").exists()

    def test_app_wraps_routes_in_boundary(self):
        src = _read("App.tsx")
        assert "ErrorBoundary" in src
        # Wrap routes, sidebar, and drawer
        assert src.count("</ErrorBoundary>") >= 3
