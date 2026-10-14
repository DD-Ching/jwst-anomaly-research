"""Offline tests for scripts/x3_g2_ident.py (X3 G2 prerequisite: astrometric identifiability)."""

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


gi = _load("x3_g2_ident")


def test_negative_mass_images_reproduce_w3_magnification():
    from jwst_anomaly import exotic_sim as es

    tau = np.array([-3.0, -2.5, -2.0, -1.0, 0.0, 1.5, 2.5, 3.0])
    A, _, _ = gi.lensed(tau, np.full_like(tau, 0.5), 0.01, n=2000)
    ref = es.light_curve(tau * 30.0, 0.0, 30.0, 0.5, 1.0, -1, 0.01)
    assert np.allclose(A, ref, rtol=5e-3, atol=1e-6)
    assert np.all(A[3:5] == 0.0)  # u < 2: no image (the W3 dip)


def test_image_centroid_points_toward_lens_and_falls_as_one_over_u():
    u = np.array([5.0, 10.0, 20.0])
    _, dx, dy = gi.lensed(u, np.zeros_like(u), 1e-4, n=50)
    assert np.all(dx < 0) and np.allclose(dy, 0.0, atol=1e-6)
    assert np.allclose(-dx * u, 1.0, rtol=0.1)


def test_blend_ratio_shift_is_affine_in_inverse_flux():
    # term 1 of any one-variable blend: x = d F_b / F, exactly a constant plus a 1/F column
    rng = np.random.default_rng(0)
    F = rng.uniform(0.3, 2.0, 50)
    x = 150.0 * 0.4 / F
    B = np.vstack([np.ones_like(F), 1.0 / F]).T
    r = x - B @ np.linalg.lstsq(B, x, rcond=None)[0]
    assert np.max(np.abs(r)) < 1e-9
