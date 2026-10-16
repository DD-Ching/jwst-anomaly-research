"""Tests for scripts/arc_diff_extended.py (compactness and local-annulus z)."""

import importlib.util
from pathlib import Path

import numpy as np

_p = Path(__file__).resolve().parents[1] / "scripts" / "arc_diff_extended.py"
_spec = importlib.util.spec_from_file_location("arc_diff_extended", _p)
ade = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ade)


def _gauss(shape, y0, x0, sig, amp):
    yy, xx = np.mgrid[: shape[0], : shape[1]]
    return amp * np.exp(-((yy - y0) ** 2 + (xx - x0) ** 2) / (2 * sig**2))


def test_compact_vs_extended_and_annulus_z():
    rng = np.random.default_rng(1)
    shape = (120, 120)
    # Smooth extended light (wide galaxy) plus one compact source elsewhere.
    e1 = _gauss(shape, 60, 60, 25.0, 1.0) + _gauss(shape, 20, 100, 1.2, 5.0)
    s = rng.normal(0.0, 1.0, shape)
    s[70, 50] = 6.0  # point brightening on the smooth light
    s[20, 100] = 6.0  # same S on the compact source
    rows = ade.measure(s, e1, np.array([[70, 50], [20, 100]], float), scale_arcsec=0.063)
    ext, comp = rows
    assert not ext["compact"] and comp["compact"]
    assert ext["s_hit"] >= 6.0
    assert 4.5 < ext["z_ann"] < 8.0
    assert ade.ext_pass(ext)
    # A wide annulus scatter (on-arc systematics) suppresses the same S.
    s2 = s.copy()
    s2[55:85, 35:65] *= 4.0
    s2[70, 50] = 6.0
    (row,) = ade.measure(s2, e1, np.array([[70, 50]], float), scale_arcsec=0.063)
    assert not ade.ext_pass(row)


def test_other_hits_are_masked_and_too_few_pixels_is_no_test():
    shape = (60, 60)
    e1 = np.ones(shape)
    s = np.random.default_rng(2).normal(0.0, 1.0, shape)
    s[30, 30] = 5.0
    s[30, 38] = 50.0  # a neighbouring hit inside the annulus
    rows = ade.measure(s, e1, np.array([[30, 30], [30, 38]], float), scale_arcsec=0.063)
    assert rows[0]["n_ann"] > 0
    assert rows[0]["z_ann"] > 3.0
    s_nan = s.copy()
    s_nan[:6, :6] = np.nan
    (edge,) = ade.measure(s_nan, e1, np.array([[1, 1]], float), scale_arcsec=0.063)
    assert not np.isfinite(edge["s_hit"]) and not ade.ext_pass(edge)
    (few,) = ade.measure(s, e1, np.array([[30, 30]], float), scale_arcsec=0.063, min_ann_px=10**6)
    assert not np.isfinite(few["z_ann"]) and ade.ext_pass(few)
