"""Unit tests for hepn_structural_stitch."""
import pytest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
from utils.hepn_structural_stitch import stitch_hepn_into_binder


def test_stitch_basic():
    """Stitch inserts HEPN1 and HEPN2 at correct positions."""
    coords = {
        "rec_end": 350,
        "hepn1_start": 360,
        "hepn1_end": 400,
        "hepn2_start": 420,
        "hepn2_end": 480,
        "linker1_len": 9,
        "linker2_len": 19,
    }
    rec = "A" * 350
    hepn1 = "RXXXH" + "Y" * 31 + "RZZZH"  # 41 aa (coords 360-400)
    hepn2 = "RAAAABBBBH" + "Z" * 51  # 61 aa (coords 420-480)
    full = rec + "X" * 9 + hepn1 + "Y" * 19 + hepn2  # 350+9+41+19+61=480
    binder = "L" * 9 + "K" * 19
    stitched = stitch_hepn_into_binder(binder, full, coords)
    assert stitched is not None
    assert stitched[:350] == rec
    assert "RXXXH" in stitched
    assert "RZZZH" in stitched
    assert len(stitched) == 480


def test_stitch_returns_none_when_baseline_too_short():
    """Returns None when baseline shorter than hepn2_end."""
    coords = {"hepn2_end": 500, "rec_end": 350, "hepn1_start": 360, "hepn1_end": 400,
              "hepn2_start": 420, "linker1_len": 9, "linker2_len": 19}
    short_full = "A" * 400
    binder = "L" * 28
    assert stitch_hepn_into_binder(binder, short_full, coords) is None
