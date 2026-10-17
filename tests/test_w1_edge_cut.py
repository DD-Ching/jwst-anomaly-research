from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("w1_edge_cut", _DIR / "w1_edge_cut.py")
wec = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(wec)


def _lims(a, b, c):
    return {m: {"upper_limit_deg2": v} for m, v in zip(wec.LIMIT_MASSES, (a, b, c), strict=True)}


def test_adopt_needs_every_limit_mass_lower():
    off = _lims(5800.0, 1100.0, 570.0)
    assert wec.adopt(off, _lims(5000.0, 1000.0, 560.0))
    assert not wec.adopt(off, _lims(5000.0, 1100.0, 560.0))  # equal is not lower
    assert not wec.adopt(off, _lims(None, 1000.0, 560.0))  # a lost limit is not lower


def test_eb_pass_needs_both_conditions():
    assert wec.eb_pass({"p_random_max": 0.01, "s_max": 4.0, "b_max": 3.5})
    assert not wec.eb_pass({"p_random_max": 0.01, "s_max": 3.4, "b_max": 3.5})
    assert not wec.eb_pass({"p_random_max": 0.07, "s_max": 4.0, "b_max": 3.5})
