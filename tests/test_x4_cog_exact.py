"""Offline tests for scripts/x4_cog_exact.py (X4 cogwheel exact-score round): the frozen fixed
points are legal and deterministic, and the script has no search command. No lalsuite or cogwheel
needed."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

_S = Path(__file__).resolve().parents[1] / "scripts"


def _load(name):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _S / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


X = _load("x4_cog_exact")
e = X.e


def test_fixed_points_legal_and_deterministic():
    pts = X.points()
    assert len(pts) == 16
    lo = np.array([bd[0] for bd in e.BOUNDS])
    hi = np.array([bd[1] for bd in e.BOUNDS])
    for x in pts.values():
        assert np.all(x >= lo) and np.all(x <= hi)
        assert x[1] <= 15 and x[2] <= 0.99 and x[5] <= 0.99
    again = X.points()
    assert all(np.array_equal(pts[k_], again[k_]) for k_ in pts)


def test_no_search_command(capsys):
    for cmd in (["refine"], ["runs"], ["one", "C1", "71"], []):
        assert X.main(cmd) == 0  # prints the usage; only `checks` computes anything
    assert "No search" in capsys.readouterr().out
