"""Unit tests for hepn_structural_stitch.

All coordinates are 0-based Python-slice convention:
    seq[hepn1_start : hepn1_end]  →  HEPN1 domain
    linker1_len = hepn1_start - rec_end
    linker2_len = hepn2_start - hepn1_end

The C-terminal tail (residues after hepn2_end) is always preserved.
"""
import pytest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
from utils.hepn_structural_stitch import stitch_hepn_into_binder


def test_stitch_basic():
    """Stitch inserts HEPN1 and HEPN2 at correct positions (0-based coords)."""
    # 0-based: hepn1 = seq[360:400] (40 aa), hepn2 = seq[420:480] (60 aa)
    coords = {
        "rec_end": 350,
        "hepn1_start": 360,
        "hepn1_end": 400,
        "hepn2_start": 420,
        "hepn2_end": 480,
        "linker1_len": 10,   # 360 - 350
        "linker2_len": 20,   # 420 - 400
    }
    rec = "A" * 350
    hepn1 = "RXXXH" + "Y" * 30 + "RZZZH"  # 40 aa (indices 360-399)
    hepn2 = "RAAAABBBBH" + "Z" * 50        # 60 aa (indices 420-479)
    c_tail = "C" * 20                        # 20 aa C-terminal tail
    full = rec + "X" * 10 + hepn1 + "Y" * 20 + hepn2 + c_tail  # 350+10+40+20+60+20=500
    binder = "L" * 10 + "K" * 20
    stitched = stitch_hepn_into_binder(binder, full, coords)
    assert stitched is not None
    assert stitched[:350] == rec
    assert "RXXXH" in stitched
    assert "RZZZH" in stitched
    assert len(stitched) == 500  # Full length including C-terminal tail
    # Verify linker regions are from binder, not baseline
    assert stitched[350:360] == "L" * 10
    assert stitched[400:420] == "K" * 20
    # Verify C-terminal tail is preserved
    assert stitched[-20:] == "C" * 20


def test_stitch_preserves_c_terminal_tail():
    """C-terminal tail after HEPN2 must always be preserved from baseline."""
    coords = {
        "rec_end": 10,
        "hepn1_start": 15,
        "hepn1_end": 25,
        "hepn2_start": 30,
        "hepn2_end": 40,
        "linker1_len": 5,   # 15 - 10
        "linker2_len": 5,   # 30 - 25
    }
    rec = "A" * 10
    hepn1 = "H" * 10         # indices 15-24
    hepn2 = "H" * 10         # indices 30-39
    c_tail = "CTERMINAL"      # 9 aa C-terminal tail
    full = rec + "X" * 5 + hepn1 + "Y" * 5 + hepn2 + c_tail  # 10+5+10+5+10+9=49
    binder = "L" * 5 + "K" * 5
    stitched = stitch_hepn_into_binder(binder, full, coords)
    assert stitched is not None
    assert stitched.endswith("CTERMINAL")
    assert len(stitched) == 49  # rec(10) + L1(5) + hepn1(10) + L2(5) + hepn2(10) + tail(9)


def test_stitch_empty_c_terminal():
    """If baseline has no residues after hepn2_end, output should match original length."""
    coords = {
        "rec_end": 350,
        "hepn1_start": 360,
        "hepn1_end": 400,
        "hepn2_start": 420,
        "hepn2_end": 480,
        "linker1_len": 10,
        "linker2_len": 20,
    }
    rec = "A" * 350
    hepn1 = "H" * 40
    hepn2 = "H" * 60
    # No C-terminal tail — baseline ends exactly at hepn2_end
    full = rec + "X" * 10 + hepn1 + "Y" * 20 + hepn2  # 350+10+40+20+60=480
    binder = "L" * 10 + "K" * 20
    stitched = stitch_hepn_into_binder(binder, full, coords)
    assert stitched is not None
    assert len(stitched) == 480  # No tail, same length as before


def test_stitch_returns_none_when_baseline_too_short():
    """Returns None when baseline shorter than hepn2_end."""
    coords = {"hepn2_end": 500, "rec_end": 350, "hepn1_start": 360, "hepn1_end": 400,
              "hepn2_start": 420, "linker1_len": 10, "linker2_len": 20}
    short_full = "A" * 400
    binder = "L" * 30
    assert stitch_hepn_into_binder(binder, short_full, coords) is None
