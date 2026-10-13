"""Tests for scripts/arcs_curvature.py (synthetic images; no network)."""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import numpy as np

_DIR = Path(__file__).resolve().parents[1] / "scripts"
_spec = importlib.util.spec_from_file_location("arcs_curvature", _DIR / "arcs_curvature.py")
ac = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ac)


def _arc_points(center_e: float, center_n: float, radius: float, n: int = 400):
    """Points on a 40-degree arc of a circle around (center_e, center_n), through the origin."""
    phi0 = np.arctan2(-center_n, -center_e)
    phi = phi0 + np.linspace(-0.35, 0.35, n)
    return center_e + radius * np.cos(phi), center_n + radius * np.sin(phi)


def test_curvature_side_points_to_the_centre():
    rng = np.random.default_rng(1)
    for ce, cn in ((5.0, 0.0), (0.0, -5.0), (-3.0, 4.0)):
        e, n = _arc_points(ce, cn, 5.0)
        e = e + rng.normal(0, 0.01, e.size)
        n = n + rng.normal(0, 0.01, n.size)
        fit = ac.fit_curvature(e, n, np.ones_like(e))
        to_centre = np.array([ce, cn]) / np.hypot(ce, cn)
        assert fit["side_e"] * to_centre[0] + fit["side_n"] * to_centre[1] > 0.95
        assert abs(abs(fit["c"]) - 0.5 / 5.0) < 0.02  # |c| = 1 / (2 R) at the vertex
        assert fit["c_err"] < 0.01


def test_classify_flat_agree_reverse():
    args = argparse.Namespace(min_sig=3.0, min_length=0.5)
    fit = {"c": 0.1, "c_err": 0.01, "side_e": 1.0, "side_n": 0.0, "length": 2.0}
    assert ac.classify(fit, 1.0, 0.0, args) == "agree"
    assert ac.classify(fit, -1.0, 0.0, args) == "reverse"
    assert ac.classify(fit | {"c_err": 0.05}, 1.0, 0.0, args) == "flat"
    assert ac.classify(fit | {"length": 0.3}, 1.0, 0.0, args) == "flat"
    assert ac.classify({}, 1.0, 0.0, args) == "no_arc"


def test_arc_pixels_ignores_halo_and_far_neighbours():
    yy, xx = np.indices((101, 101))
    halo = 50.0 * np.exp(-((xx - 100) ** 2 + (yy - 100) ** 2) / (2 * 60.0**2))
    arc = 20.0 * np.exp(-((yy - 50) ** 2) / 4.0) * (np.abs(xx - 50) < 15)
    neighbour = 30.0 * np.exp(-((xx - 90) ** 2 + (yy - 50) ** 2) / 8.0)
    img = halo + arc + neighbour + np.random.default_rng(2).normal(0, 0.5, halo.shape)
    mask, _ = ac.arc_pixels(img, 3.0, 1.0, radius_px=25.0, bkg_px=31)
    assert mask[50, 50] and mask[50, 60]
    assert not mask[50, 90]  # outside the window
    assert not mask[10, 10]  # halo removed


def test_image_uri():
    assert ac.image_uri("jw01208-o002_t001_nircam_clear-f200w_i2d.fits") == (
        "s3://stpubdata/jwst/public/jw01208/L3/t/o002/jw01208-o002_t001_nircam_clear-f200w_i2d.fits"
    )
