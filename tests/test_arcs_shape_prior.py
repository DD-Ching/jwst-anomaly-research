"""Rule functions of scripts/arcs_shape_prior.py (written before the run)."""

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("arcs_shape_prior", _DIR / "arcs_shape_prior.py")
sp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sp)


def test_eps_from_e():
    assert sp.eps_from_e(np.array([0.0, 0.5, 1.0])) == pytest.approx([0.0, 1 / 3, 1.0])


def test_prior_cuts():
    t = Table(
        {
            "is_extended": [True, True, False, True, True],
            "semimajor_px": [3.0, 1.0, 3.0, 3.0, 3.0],
            "snr": [20.0, 20.0, 20.0, 5.0, 20.0],
            "ellipticity": [0.5, 0.5, 0.5, 0.5, np.nan],
        }
    )
    assert sp.prior_eps(t, 2.0, 10.0) == pytest.approx([1 / 3])


def test_p_anti_isotropic_at_zero_shear():
    eps = np.random.default_rng(0).uniform(0.34, 0.9, 5000)
    assert sp.p_anti_emp(0.0, eps, n=100_000) == pytest.approx(1 / 3, abs=0.01)


def test_p_anti_matches_gaussian_code():
    """With a Gaussian-drawn prior the empirical code reproduces `anti_probability`."""
    rng = np.random.default_rng(3)
    eps = np.hypot(rng.normal(0, 0.26, 400_000), rng.normal(0, 0.26, 400_000))
    eps = eps[eps < 0.95]
    for g in (0.2, 0.5, 1.5):
        assert sp.p_anti_emp(g, eps) == pytest.approx(sp.lc.anti_probability(g), abs=0.01)


def test_p_anti_falls_with_shear():
    eps = np.random.default_rng(1).uniform(0.0, 0.9, 5000)
    assert sp.p_anti_emp(0.6, eps) < sp.p_anti_emp(0.2, eps) < 1 / 3


def test_gate_and_outcome():
    assert sp.gate(9, 9.0, 1 / 3)["passes"]
    assert not sp.gate(9, 9.0, 0.36)["passes"]
    assert not sp.gate(20, 5.0, 1 / 3)["passes"]
    assert sp.outcome(False, 0.9) == "tool fails; undetermined"
    assert sp.outcome(True, 0.05) == "pooled excess explained by the empirical shape prior"
    assert sp.outcome(True, 0.049) == "excess persists under the empirical prior"
