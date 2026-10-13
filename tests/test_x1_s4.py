"""Offline tests for scripts/x1_s4.py (model-free far-side disk bound), on synthetic data only."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import x1_s4  # noqa: E402


def _data(xy, sigma=10.0, flux=None):
    xy = np.asarray(xy, float)
    flux = np.full(len(xy), 0.5) if flux is None else np.asarray(flux, float)
    return {"xy": xy, "err": np.full_like(xy, sigma), "flux": flux}


def test_points_inside_disk_fit_exactly():
    ang = np.linspace(0, 2 * math.pi, 10, endpoint=False)
    xy = 15.0 * np.c_[np.cos(ang), np.sin(ang)] + [100.0, -40.0]
    b = x1_s4.disk_bound(_data(xy), r=20.0, d_cap=20.0, fs_max=0.0)
    assert b["chi2_lb"] < 1e-6


def test_two_opposite_points_analytic():
    # separation 2 (r + delta), isotropic sigma, no static light: chi2 = 2 delta^2 / sigma^2
    r, delta, sigma = 20.0, 7.0, 10.0
    xy = [[r + delta, 0.0], [-(r + delta), 0.0]]
    b = x1_s4.disk_bound(_data(xy, sigma), r=r, d_cap=20.0, fs_max=0.0)
    assert b["chi2_lb"] == pytest.approx(2 * delta**2 / sigma**2, rel=1e-5)
    assert b["epochs_on_disk_edge"] == 2


def test_bound_never_exceeds_a_track_inside_the_disk():
    # any concrete far-side-like track (inside the disk) plus static light fits at least as well
    rng = np.random.default_rng(1)
    t = np.linspace(0, 1, 10)
    flux = 0.4 + 0.2 * t
    u = 18.0 * np.c_[np.cos(5 * t), np.sin(5 * t)]
    fs, off = 0.3, np.array([15.0, -5.0])
    s = fs * flux.mean()
    xy = ((flux - s)[:, None] * u + s * off) / flux[:, None] + rng.normal(0, 10, (10, 2))
    d = _data(xy, flux=flux)
    w = 1.0 / d["err"] ** 2
    model = ((flux - s)[:, None] * u + s * off) / flux[:, None]
    res = xy - model
    c0 = (res * w).sum(0) / w.sum(0)
    chi2_true = float((((res - c0) ** 2) * w).sum())
    b = x1_s4.disk_bound(d, r=20.0, d_cap=20.0)
    assert b["chi2_lb"] <= chi2_true + 1e-6


def test_static_light_offset_is_used():
    # a large off-centre static source with a 1/F_t signature cannot be fit with fs_max = 0
    flux = np.linspace(0.3, 1.0, 10)
    q = np.array([60.0, 0.0])
    xy = q / flux[:, None]
    loose = x1_s4.disk_bound(_data(xy, flux=flux), r=1.0, d_cap=400.0)
    tight = x1_s4.disk_bound(_data(xy, flux=flux), r=1.0, d_cap=400.0, fs_max=0.0)
    assert loose["chi2_lb"] < 1e-4 < tight["chi2_lb"]
