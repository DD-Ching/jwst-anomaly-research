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


def test_s1c_fast_fft_length_never_coarser():
    e = _load("x4_s1c_r1")
    for pad in (32, 64):
        assert e._nfft(pad) >= pad * e.b._FREQS.size


def test_s1c_grid_scan_matches_pointwise_overlap(monkeypatch):
    e = _load("x4_s1c_r1")
    n = e.b._FREQS.size
    monkeypatch.setattr(e.b, "weight", lambda: e.b._BAND.astype(float))
    rng = np.random.default_rng(1)
    f = np.exp(1j * rng.uniform(0, 1, (3, n))) * rng.uniform(0.5, 2, (3, n))
    monkeypatch.setitem(e._G, "F", f)
    monkeypatch.setitem(e._G, "F2", np.abs(f) ** 2)
    monkeypatch.setitem(e._G, "nodes", [(1.0, 0.1), (2.0, 0.2), (3.0, 0.3)])
    t = (rng.normal(size=n) + 1j * rng.normal(size=n)) * e.b._BAND
    h = (rng.normal(size=n) + 1j * rng.normal(size=n)) * e.b._BAND
    for fam in ("A", "D"):
        m = e.grid_scan(t, h, fam)
        for k in range(3):
            assert abs(m[k] - e.ov(t, h * f[k], fam)) < 1e-12


def test_s1c_domain_and_bounds():
    e = _load("x4_s1c_r1")
    lmc = np.log(e.MC0)
    assert e.in_domain([lmc, 3.0, 0.5, -0.5], "A")
    assert not e.in_domain([lmc, 16.0, 0.5, -0.5], "A")
    assert not e.in_domain([lmc, 3.0, 0.8, 0.0, 0.7, 0, 0, 0, 1.0, 0, 0], "D")
    assert e.bounds([lmc, 15.0, 0.99, 0.0], "A", (2441.0, 0.096)) == ["q_max", "chi"]


def test_s1d_polarisation_max_equals_psi_scan(monkeypatch):
    e = _load("x4_s1d_de")
    n = e.b._FREQS.size
    monkeypatch.setattr(e.b, "weight", lambda: e.b._BAND.astype(float))
    rng = np.random.default_rng(5)
    hp, hc, t = ((rng.normal(size=n) + 1j * rng.normal(size=n)) * e.b._BAND for _ in range(3))
    t = t + 0.8 * (np.cos(0.6) * hp + np.sin(0.6) * hc)
    scan = max(
        e.m.ov(t, np.cos(2 * p) * hp + np.sin(2 * p) * hc, "D") for p in np.linspace(0, np.pi, 361)
    )
    val = e.stat(t, hp, hc)
    assert val >= scan - 1e-12 and val - scan < 1e-3


def test_s1d_spin_parametrisation_round_trip():
    e = _load("x4_s1d_de")
    x = np.array([np.log(e.MC0), 2.0, 0.5, 0.3, 1.0, 0.7, -0.4, 4.0, 1.1, 2.0])
    assert np.allclose(e.from_cart(*e.to_cart(x)), x, atol=1e-12)
    for lo_hi, v in zip(e.BOUNDS, x, strict=True):
        assert lo_hi[0] <= v <= lo_hi[1]
    assert e.bounds_hit(x) == []
