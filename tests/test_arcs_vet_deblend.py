"""Tests for scripts/arcs_vet_deblend.py (synthetic values; no network)."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np

_DIR = Path(__file__).resolve().parents[1] / "scripts"
_spec = importlib.util.spec_from_file_location("arcs_vet_deblend", _DIR / "arcs_vet_deblend.py")
vd = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(vd)


def test_deblend_class():
    assert vd.deblend_class(10, 3, 3, 0.10) == "two-knot"
    assert vd.deblend_class(10, 2.9, 3, 0.10) == "undecided"
    assert vd.deblend_class(10, 3, 3, 0.09) == "undecided"
    assert vd.deblend_class(-10, 9, 9, 0.3) == "single"
    assert vd.deblend_class(float("nan"), 9, 9, 0.3) == "undecided"


def test_gate_ok():
    assert vd.gate_ok(0.8, 0.8)
    assert not vd.gate_ok(0.75, 1.0)


def test_sed_match_and_floor():
    f = np.array([1.0, 2.0, 3.0, 2.0])
    e = np.full(4, 0.01)
    p, n = vd.sed_match(f, e, 0.5 * f, e)
    assert n == 4 and p > 0.5
    g = f * np.array([1.0, 1.03, 0.97, 1.0])
    assert vd.sed_class(vd.sed_match(f, e, g, e)[0]) == "SED-different"
    assert vd.sed_class(vd.sed_match(f, e, g, e, floor=0.1)[0]) == "SED-consistent"
    assert np.isnan(vd.sed_match(f[:1], e[:1], f[:1], e[:1])[0])


def test_ratio_consistency_drops_faint_bands():
    f1 = np.array([1.0, 1.0, 1.0])
    f2 = np.array([2.0, 2.0, 0.01])
    e = np.full(3, 0.05)
    p, n = vd.ratio_consistency(f1, e, f2, e)
    assert n == 2 and p > 0.5


def test_fit_models_recovers_two_knots():
    rng = np.random.default_rng(1)
    yy, xx = np.indices((41, 41), dtype=float)
    img = vd._gauss(xx, yy, 17, 20, 50, 1.3) + vd._gauss(xx, yy, 23, 20, 50, 1.3)
    img += rng.normal(0, 0.05, img.shape)
    f = vd.fit_models(img, np.full(img.shape, 0.05), 20, 20, 12)
    assert f["dbic"] > 10 and abs(f["sep_px"] - 6) < 0.3
    single = vd._gauss(xx, yy, 20, 20, 100, 3.5, 0.45, 0.3) + rng.normal(0, 0.05, img.shape)
    g = vd.fit_models(single, np.full(img.shape, 0.05), 20, 20, 12)
    assert g["dbic"] < 10


def test_gate_uses_requested_separation():
    rng = np.random.default_rng(2)
    px = 0.031
    sci = rng.normal(0, 0.05, (320, 320))
    err = np.full(sci.shape, 0.05)
    sizes = {"single": (40.0, 2.0, 0.6), "two": (20.0, 1.3, 20.0, 1.3)}
    g = vd.gate(sci, err, 160, 160, px, 0.4 / px, sizes, 3, rng, sep_arcsec=0.25)
    assert g["n_single"] == 3 and g["n_two"] == 3
    assert g["frac_two_called_two_knot"] == 1.0 and g["passed"]
