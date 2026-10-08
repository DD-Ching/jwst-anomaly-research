"""Offline tests for scripts/w5_counts.py helpers (D-063): null model and lens geometry."""

from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import numpy as np
import pytest

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("w5_counts", _DIR / "w5_counts.py")
w5 = importlib.util.module_from_spec(_spec)
sys.modules["w5_counts"] = w5
_spec.loader.exec_module(w5)


def test_null_model_exponential_tail():
    """An exponential peak distribution is recovered and z_det matches the requested FAP."""
    rng = np.random.default_rng(0)
    z = rng.exponential(0.5, 200_000)
    m = w5.null_model(z, n_search=200_000)
    assert m["tau"] == pytest.approx(0.5, rel=0.05)
    assert m["z_flag"] == pytest.approx(z.max())
    # expected false count: N exp(-z / tau)
    assert w5.n_false(m, 5.0) == pytest.approx(2e5 * math.exp(-10.0), rel=0.3)
    assert w5.n_false(m, m["z_det"]) <= w5.FAP * 1.0001
    assert m["z_det"] >= m["z_flag"]


def test_theta_mass_geometry():
    """θ_E = 2.85 mas (M/M☉)^½ at 1 kpc (D_S >> D_L); θ_E ∝ M^½; inverse round trip."""
    g = w5.GEOMETRIES[0]
    th = w5.theta_e_arcmin([1.0, 4.0], g)
    assert th[0] * 60e3 == pytest.approx(2.85, rel=0.01)
    assert th[1] == pytest.approx(2 * th[0])
    for geo in w5.GEOMETRIES:
        m = w5.mass_for_theta(10.0, geo)
        assert w5.theta_e_arcmin(m, geo)[0] == pytest.approx(10.0)
    # a cosmological lens needs more mass for the same θ_E than one at 1 Mpc
    assert w5.mass_for_theta(1.0, w5.GEOMETRIES[2]) > w5.mass_for_theta(1.0, w5.GEOMETRIES[1])


def test_regions_are_disjoint_and_chunked():
    a, b = w5.REGIONS.values()
    assert a.ra_max <= b.ra_min
    assert len(a.chunks()) == 50 and len(b.chunks()) == 50
