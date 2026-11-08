"""Offline tests of the GLAFIC-realisation route for RXC J0600 581 (D-024 addendum 18)."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
from astropy.wcs import WCS

_P = Path(__file__).resolve().parents[1] / "scripts" / "arcs_shear_error_581g.py"
_spec = importlib.util.spec_from_file_location("arcs_shear_error_581g", _P)
m = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m)


def _wcs(n: int, scale_arcsec: float, ra=90.0, dec=-20.0) -> WCS:
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [ra, dec]
    w.wcs.crpix = [n / 2, n / 2]
    s = scale_arcsec / 3600.0
    w.wcs.cd = [[-s, 0.0], [0.0, s]]
    return w


def test_pixel_scale():
    assert np.isclose(m.pixel_scale_arcsec(_wcs(10, 0.4)), 0.4)


def test_resample_linear_field_exact():
    fine = _wcs(400, 0.1)
    coarse = _wcs(100, 0.4)
    jj, ii = np.mgrid[0:400, 0:400]
    ra, dec = fine.pixel_to_world_values(ii, jj)
    data = (3.0 * (ra - 90.0) + 5.0 * (dec + 20.0)).astype(np.float64)
    out = m.resample(data, fine, coarse, (100, 100))
    cj, ci = np.mgrid[0:100, 0:100]
    cra, cdec = coarse.pixel_to_world_values(ci, cj)
    want = 3.0 * (cra - 90.0) + 5.0 * (cdec + 20.0)
    ok = np.isfinite(out)
    assert ok.mean() > 0.9
    assert np.allclose(out[ok], want[ok], atol=1e-6)


def test_resample_outside_is_nan():
    fine = _wcs(40, 0.1)
    coarse = _wcs(40, 0.4)
    out = m.resample(np.ones((40, 40)), fine, coarse, (40, 40))
    assert np.isnan(out[0, 0]) and np.isfinite(out[20, 20])


def test_v1_v2_gates():
    assert m.v1_gate(5e-4) and not m.v1_gate(2e-3) and not m.v1_gate(float("nan"))
    assert m.v2_gate(["aligned", "anti"], ["aligned", "anti"], [0.1, 1.9])
    assert not m.v2_gate(["aligned", "anti"], ["aligned", "anti"], [0.1, 2.1])
    assert not m.v2_gate(["aligned", "anti"], ["aligned", "mixed"], [0.1, 0.1])
    assert not m.v2_gate([], [], [])


def test_v3_stat_and_gate():
    best = np.array([[10.0, 0.0], [0.0, 20.0], [5.0, 5.0]])
    rng = np.random.default_rng(1)
    real = best[None] + rng.normal(0.0, 0.01, (99, 3, 2))
    s = m.v3_stat(best, real)
    assert s < 0.01 and m.v3_gate(s)
    # a factor-4 unit error fails
    assert not m.v3_gate(m.v3_stat(best, 4.0 * real))
    assert not m.v3_gate(float("nan"))


def test_frozen_constants():
    assert m.N_RANGE == 99 and m.RANGE_GATE_MAPS == (0, 49, 98)
    assert m.ROW == 581 and m.Z_LENS == 0.43


def test_same_grid_includes_rotation():
    ref = _wcs(10, 0.4).to_header()
    for k, v in (("NAXIS1", 10), ("NAXIS2", 10)):
        ref[k] = v
    assert m.same_grid(ref.copy(), ref)
    rot = ref.copy()
    rot["PC1_2"] = 0.0
    rot["CD1_2"] = 1e-6
    ref["CD1_2"] = 0.0
    assert not m.same_grid(rot, ref)
    shifted = ref.copy()
    shifted["CRVAL1"] = ref["CRVAL1"] + 1e-3
    assert not m.same_grid(shifted, ref)


def test_local_peak_interior_and_edge():
    w = _wcs(101, 0.1)
    jj, ii = np.mgrid[0:101, 0:101]
    k = 0.01 * ii.astype(float)  # rises to +x: box maximum on the edge
    k[52, 50] = 5.0  # interior peak 0.2 arcsec from the centre pixel (50, 50)
    ra, dec = w.pixel_to_world_values(50, 50)
    out = m.local_peak(k, w, float(ra), float(dec), 2.0)
    assert out["kappa_max_box"] == 5.0 and np.isclose(out["max_dist_arcsec"], 0.2)
    k[52, 50] = 0.0
    out = m.local_peak(k, w, float(ra), float(dec), 2.0)
    assert out["max_dist_arcsec"] >= 2.0
