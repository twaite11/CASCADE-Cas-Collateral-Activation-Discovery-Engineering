"""Canonical CRISPR-array repeat / spacer length bounds.

B-21 fix: previously these values were duplicated across
`fix_crrna_assignments.py`, `discover_crrna.py` (as magic numbers
"len >= 20", "len >= 23"), and the mining scripts' implicit ranges in
their k-mer searches.  When the audit tightened the cas13a DR floor from
20 nt to 23 nt to drop spurious tRNAs, only `fix_crrna_assignments.py`
got updated -- `discover_crrna.py:262` still accepted 20-nt "repeats".

This module is the single source of truth.  Importers:
  * scripts/fix_crrna_assignments.py
  * scripts/discover_crrna.py
  * scripts/mining_v3.py (via fix_crrna_assignments.find_crispr_arrays)

The numeric choices below match the literature consensus for Cas13a/b/c/d
direct repeats (28-37 nt for class 2 type VI; 23-50 nt umbrella range).
"""
from __future__ import annotations

# Direct-repeat length bounds (nt). The umbrella that covers every known
# class 2 type VI subtype with margin on either side.
CRISPR_REPEAT_MIN = 23
CRISPR_REPEAT_MAX = 50

# Spacer length bounds (nt).  Cas13 spacers cluster around 28-32 nt; we
# widen to 15-80 to also catch class 1 / type III arrays that may
# coincidentally hit the DR length window during pre-filter scans.
CRISPR_SPACER_MIN = 15
CRISPR_SPACER_MAX = 80

# Minimum repeats per array. 3 is the canonical "real CRISPR array"
# threshold; using 2 lets a tRNA inverted-repeat masquerade as an array.
MIN_ARRAY_UNITS = 3

__all__ = [
    "CRISPR_REPEAT_MIN",
    "CRISPR_REPEAT_MAX",
    "CRISPR_SPACER_MIN",
    "CRISPR_SPACER_MAX",
    "MIN_ARRAY_UNITS",
]
