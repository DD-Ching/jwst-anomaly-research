"""Tests for scripts/w5_euclid_psf_radial.py (synthetic inputs only)."""

from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import numpy as np

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location(
    "w5_euclid_psf_radial", _DIR / "w5_euclid_psf_radial.py"
)
pr = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(pr)


def _grid(n=60, half=300.0):
    g = np.linspace(-half, half, n)
    x, y = np.meshgrid(g, g)
    return x.ravel(), y.ravel()


def test_radial_pattern_sign_and_amplitude():
    # A pure radial pattern of amplitude 0.02 about the origin, built with w5.radial_shear's phase.
    x, y = _grid()
    phase = pr.w5.radial_shear(x, y, 0.0, 0.0, 1.0, 0.0)
    e = 0.02 * phase / np.abs(phase)
    r = pr.annulus_radial(x, y, e, [0.0], [0.0], 90.0, 180.0)
    assert math.isclose(r["mean"][0], 0.02, rel_tol=1e-9)
    t = pr.annulus_radial(x, y, -e, [0.0], [0.0], 90.0, 180.0)  # tangential is negative
    assert math.isclose(t["mean"][0], -0.02, rel_tol=1e-9)


def test_constant_anisotropy_cancels():
    x, y = _grid()
    e = np.full(len(x), 0.015 - 0.01j)
    r = pr.annulus_radial(x, y, e, [0.0], [0.0], 90.0, 180.0)
    assert abs(r["mean"][0]) < 1e-3  # a constant PSF ellipticity averages out over a full annulus


def test_min_stars_skips_sparse_centres():
    r = pr.annulus_radial([100.0, 0.0], [0.0, 120.0], [0.01, 0.01], [0.0], [0.0], 90.0, 180.0)
    assert np.isnan(r["mean"][0]) and r["n"][0] == 2


def test_signal_and_decide():
    assert math.isclose(pr.SIGNAL, 0.5 * 2 * math.log(2) / 6.75)

    def c(r1, chi2):
        return {"rho": {"1": r1, "2.3": 2.3 * r1}, "chi2_over_dof": chi2}

    assert pr.decide([c(0.05, 3.0), c(0.08, 0.5)]) == "PSF gradients cannot fake or hide a flag"
    # the largest rho sits in a cell with chi2 <= dof; another cell still meets rule 2
    assert pr.decide([c(1.2, 0.98), c(0.5, 1.01)]) == (
        "PSF term not negligible at the edge of the resolved cut"
    )
    assert pr.decide([c(1.0, 1.5), c(0.5, 1.5)]) == "PSF could fake or hide a flag"
    assert pr.decide([c(1.2, 0.9), c(0.5, 0.95)]).startswith("no PSF pattern above noise")


def test_error_uses_field_rms():
    x, y = _grid()
    e = np.where(np.arange(len(x)) % 2, 0.03 + 0j, -0.03 + 0j)
    r = pr.annulus_radial(x, y, e, [0.0], [0.0], 90.0, 180.0)
    assert math.isclose(r["err"][0], math.sqrt(0.03**2 / 2) / math.sqrt(r["n"][0]))


def test_annulus_counts():
    x, y = np.array([0.0, 100.0, 150.0, 200.0]), np.zeros(4)
    n = np.array([5.0, 2.0, 3.0, 7.0])
    np.testing.assert_allclose(pr.annulus_counts(x, y, n, [0.0], [0.0], 90.0, 180.0), [5.0])


def test_coverage_query_groups_cells():
    q = pr.coverage_query("EDF-F")
    assert "GROUP BY iy, ix" in q and "point_like_prob < 0.1" in q and "FLOOR(declination*120)" in q
