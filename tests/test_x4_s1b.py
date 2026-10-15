"""Offline tests for scripts/x4_s1b_checks.py and x4_s1b_xphm.py (X4 S1b). No lalsuite needed."""

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


c = _load("x4_s1b_checks")
d = _load("x4_s1b_xphm")


def test_nm_restart_converges_on_quadratic():
    r = c.nm_restart(lambda x: float(np.sum((x - 1.0) ** 2)), [0.0, 0.0], np.array([0.5, 0.5]))
    assert r["converged_restart"] and r["f"] < 1e-8
    assert np.allclose(r["x"], 1.0, atol=1e-3)
    assert r["history"] == sorted(r["history"], reverse=True)


def test_bounds_hit_flags():
    assert c.bounds_hit([0.0, 15.0, 0.2, 0.2], 15.0) == ["q_max"]
    assert c.bounds_hit([0.0, 3.0, 0.99, 0.0], 15.0) == ["chi"]
    assert c.bounds_hit([0.0, 3.0, 0.1, 0.1, np.log(500.0), 1.0], 15.0) == []


def test_frozen_representative_cases():
    assert c.CASES["R1"] == (30.0, 300.0, 0.3) and c.CASES["R4"] == (200.0, 100.0, 1.0)
    assert d.D_CASES == ("R2", "R3", "R4")


def test_overlap_real_has_no_free_phase(monkeypatch):
    monkeypatch.setattr(d.b, "weight", lambda: np.ones(d.b._FREQS.size))
    rng = np.random.default_rng(0)
    n = c.b._FREQS.size
    h = (rng.normal(size=n) + 1j * rng.normal(size=n)) * c.b._BAND
    assert abs(d.overlap_real(h, 2.5 * h) - 1.0) < 1e-9
    # a 90 degree rotation is only partly absorbed (through the time shift), unlike |corr|
    assert d.overlap_real(h * 1j, h) < 0.9 < c.b.overlap(h * 1j, h)
