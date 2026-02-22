"""
HEPN sequence stitching for PXDesign output.
PXDesign generates linkers only (REC→HEPN1 and HEPN1→HEPN2). This module inserts
the wild-type HEPN1 and HEPN2 sequences at the correct indices into the designed binder.
"""
import logging

log = logging.getLogger(__name__)


def stitch_hepn_into_binder(
    binder_seq: str,
    baseline_full_seq: str,
    coords: dict,
) -> str:
    """
    Insert wild-type HEPN1 and HEPN2 into the designed binder sequence.
    The binder has no native gaps; we precisely index where to insert the HEPN domains.

    full = rec + binder[:linker1_len] + hepn1 + binder[linker1_len:linker1_len+linker2_len] + hepn2

    Returns stitched full sequence (rec + linkers + hepn1 + linkers + hepn2), or None on failure.
    """
    rec_end = coords["rec_end"]
    linker1_len = coords["linker1_len"]
    linker2_len = coords["linker2_len"]
    hepn1_start = coords["hepn1_start"]
    hepn1_end = coords["hepn1_end"]
    hepn2_start = coords["hepn2_start"]
    hepn2_end = coords["hepn2_end"]

    expected_binder_len = linker1_len + linker2_len
    if len(baseline_full_seq) < hepn2_end:
        log.warning(f"Baseline sequence too short for HEPN2 end ({hepn2_end})")
        return None

    hepn1_seq = baseline_full_seq[hepn1_start - 1 : hepn1_end]
    hepn2_seq = baseline_full_seq[hepn2_start - 1 : hepn2_end]
    rec_seq = baseline_full_seq[:rec_end]

    if expected_binder_len == 0:
        # No native linkers; treat entire binder as connector between REC and HEPN1
        linker1 = binder_seq
        linker2 = ""
    elif len(binder_seq) < expected_binder_len:
        linker1 = binder_seq[: min(linker1_len, len(binder_seq))]
        remainder = binder_seq[len(linker1) :]
        linker2 = remainder[: linker2_len] if len(remainder) >= linker2_len else remainder + "G" * max(0, linker2_len - len(remainder))
        log.warning(f"Binder shorter than expected ({len(binder_seq)} < {expected_binder_len}); using available residues")
    elif len(binder_seq) > expected_binder_len:
        linker1 = binder_seq[:linker1_len]
        linker2 = binder_seq[linker1_len : linker1_len + linker2_len]
    else:
        linker1 = binder_seq[:linker1_len]
        linker2 = binder_seq[linker1_len : linker1_len + linker2_len]

    full = rec_seq + linker1 + hepn1_seq + linker2 + hepn2_seq
    return full
