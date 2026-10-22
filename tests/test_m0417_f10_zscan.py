"""Frozen rules of ``scripts/m0417_f10_zscan.py`` (MACS0417 family-10 z-scan)."""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import m0417_f10_zscan as fz  # noqa: E402


def _row(z, ok, sep3):
    return {"z": z, "lands": ok, "sep_third_arcsec": sep3}


def test_lands_needs_three_finite_within_match():
    assert fz.lands([0.2, 0.3, 1.0])
    assert not fz.lands([0.2, 0.3, 1.01])
    assert not fz.lands([0.2, 0.3])
    assert not fz.lands([0.2, math.nan, 0.1])
    assert not fz.lands([])


def test_known_ok_pair():
    assert fz.known_ok([0.1, 0.9])
    assert not fz.known_ok([0.1, 1.2])
    assert not fz.known_ok([math.nan, math.nan])


def test_outcome_classes():
    assert fz.outcome(False, [_row(1.8, True, 0.1)])["outcome"].startswith("not read")
    a = fz.outcome(True, [_row(1.7, False, 1.4), _row(1.8, True, 0.5), _row(2.5, True, 0.6)])
    assert a["outcome"].startswith("A:")
    assert a["z_lands"] == [1.8, 2.5] and a["z_lands_in_photoz"] == [1.8, 1.8]
    a2 = fz.outcome(True, [_row(1.8, False, 1.4), _row(2.5, True, 0.6)])
    assert a2["outcome"].startswith("A'") and a2["z_lands_in_photoz"] is None
    b = fz.outcome(True, [_row(1.8, False, 1.4), _row(2.5, False, None), _row(3.0, False, 1.1)])
    assert b["outcome"].startswith("B") and b["min_sep_third"] == {"z": 3.0, "sep_arcsec": 1.1}


def test_photoz_edges_inclusive():
    assert fz.outcome(True, [_row(1.72, True, 0.2)])["outcome"].startswith("A:")
    assert fz.outcome(True, [_row(1.88, True, 0.2)])["outcome"].startswith("A:")
