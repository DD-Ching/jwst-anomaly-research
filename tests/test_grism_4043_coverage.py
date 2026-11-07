"""Offline tests of the footprint geometry in scripts/grism_4043_coverage.py."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import grism_4043_coverage as g  # noqa: E402

SQUARE = "POLYGON 110.0 -73.0 110.1 -73.0 110.1 -73.1 110.0 -73.1"


def test_parse_two_polygons():
    polys = g.parse_polygons(SQUARE + " " + SQUARE.replace("110.0", "111.0"))
    assert len(polys) == 2
    assert polys[0].shape == (4, 2)


def test_inside_and_outside_high_dec():
    (p,) = g.parse_polygons(SQUARE)
    assert g.inside(110.05, -73.05, p)
    assert not g.inside(110.15, -73.05, p)
    assert not g.inside(110.05, -72.99, p)


def test_edge_distance():
    (p,) = g.parse_polygons(SQUARE)
    # 0.01 deg in Dec from the northern edge = 36 arcsec.
    assert np.isclose(g.edge_distance_arcsec(110.05, -73.01, p), 36.0, rtol=1e-3)


def test_coverage_reports_nearest_edge_only_when_inside():
    fp = {"a": g.parse_polygons(SQUARE)}
    assert g.coverage(110.05, -73.01, fp)["a"]["inside"]
    assert g.coverage(112.0, -73.01, fp)["a"] == {"inside": False, "edge_arcsec": None}


import grism_4043_lines as gl  # noqa: E402


def _spec(seed=0):
    rng = np.random.default_rng(seed)
    wave = np.linspace(3.9, 5.0, 1100)
    return wave, 1e-6 + rng.normal(0, 1e-8, wave.size)


def test_line_snr_finds_injected_line_and_not_noise():
    wave, flux = _spec()
    c = gl.LINES_UM["Pabeta"] * 3.84
    s, w = gl.line_snr(wave, gl.inject(wave, flux, c, 10.0), gl.LINES_UM["Pabeta"])
    assert s > 7 and abs(w - c) < 0.002
    s0, _ = gl.line_snr(wave, flux, gl.LINES_UM["Pabeta"])
    assert s0 < 5


def test_decide_needs_both_grisms_and_matching_wavelength():
    lines = {ln: (1.0, 4.5) for ln in gl.LINES_UM}
    both = {"r": dict(lines, Pabeta=(6.0, 4.92)), "c": dict(lines, Pabeta=(6.0, 4.921))}
    assert gl.decide(both)["verdict"] == "emission_detected"
    one = {"r": dict(lines, Pabeta=(6.0, 4.92)), "c": dict(lines, Pabeta=(4.0, 4.92))}
    assert gl.decide(one)["verdict"] == "no_line_detected"
    apart = {"r": dict(lines, Pabeta=(6.0, 4.92)), "c": dict(lines, Pabeta=(6.0, 4.93))}
    assert not gl.decide(apart)["detected"]["Pabeta"]
    nan = {
        "r": dict(lines, Pabeta=(float("nan"), float("nan"))),
        "c": dict(lines, Pabeta=(6.0, 4.92)),
    }
    assert not gl.decide(nan)["detected"]["Pabeta"]


def test_inject_off_spectrum_leaves_flux_finite_and_unchanged():
    wave, flux = _spec()
    out = gl.inject(wave, flux, 5.3, 7.0)
    assert np.array_equal(out, flux)


def test_run_injections_counts_coverage():
    wave, flux = _spec()
    short = (wave[wave < 4.76], flux[wave < 4.76])
    res = gl.run_injections({"r": (wave, flux), "c": short}, n=20, snr=12.0, seed=0)
    assert res["Pabeta"]["n_both_covered"] < 20
    assert np.isfinite(res["Pabeta"]["recovery"])
