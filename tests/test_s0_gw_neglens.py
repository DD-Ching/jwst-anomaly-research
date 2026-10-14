"""Offline tests for scripts/s0_gw_neglens.py (S0-2026-10-15 A: GW lensing by a negative mass)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
from scipy.special import jv

_P = Path(__file__).resolve().parents[1] / "scripts" / "s0_gw_neglens.py"
_spec = importlib.util.spec_from_file_location("s0_gw_neglens", _P)
gl = importlib.util.module_from_spec(_spec)
sys.modules.setdefault("s0_gw_neglens", gl)
_spec.loader.exec_module(gl)


def _numeric(w, y, s, n=200000):
    """Diffraction integral with the contour rotated by pi/4 (converges for both signs)."""
    r = np.linspace(1e-9, 12 / np.sqrt(w) + 2 * y, n)
    x = r * np.exp(1j * np.pi / 4)
    f = x * jv(0, w * x * y) * np.exp(1j * w * (x * x / 2 - s * np.log(x))) * np.exp(1j * np.pi / 4)
    return -1j * w * np.exp(1j * w * y * y / 2) * np.trapezoid(f, r)


def test_closed_form_matches_diffraction_integral_both_signs():
    for s in (1, -1):
        for w, y in ((0.5, 0.3), (2.0, 0.5), (5.0, 1.0)):
            assert abs(gl.amp_factor(w, y, s) / _numeric(w, y, s) - 1) < 1e-5


def test_low_frequency_limit_and_negative_mass_shadow():
    assert abs(gl.amp_factor(1e-4, 0.5, -1) - 1) < 1e-3
    # behind a negative-mass lens (y < 2) high frequencies are suppressed: the wave-optics W3 dip
    assert abs(gl.amp_factor(10.0, 1.0, -1)) < 0.1
    assert abs(gl.amp_factor(10.0, 1.0, +1)) > 1.0


def test_geometric_limit_matches_w3_magnification():
    from jwst_anomaly import exotic_sim as es

    w = np.linspace(200.0, 260.0, 400)  # average over interference fringes
    mag = np.mean([abs(gl.amp_factor(x, 3.0, -1)) ** 2 for x in w])
    ref = float(es.light_curve(np.array([0.0]), 0.0, 1.0, 3.0, 1.0, -1, 0.0)[0])
    assert abs(mag / ref - 1) < 0.03
