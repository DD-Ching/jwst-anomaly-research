"""Offline tests for scripts/s3_hybrid.py (S3 hybrid images in lensed-quasar light curves)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("s3_hybrid", _DIR / "s3_hybrid.py")
s3 = importlib.util.module_from_spec(_spec)
sys.modules["s3_hybrid"] = s3  # dataclasses look their module up
_spec.loader.exec_module(s3)


def test_lag_window_follows_the_sis_leg_split():
    sy = s3.System("x", "lcab/x.dat", -100.0, 5.0, 0.5, 2.0)
    w = s3.lag_window(sy)
    q = w["dls_over_ds"]
    assert 0.0 < q < 1.0
    assert w["lead"] == "A" and w["trail"] == "B"
    assert np.allclose(w["into_trail"], (-100.0 * (q + 1.0), -100.0 * q))
    assert np.allclose(w["into_lead"], (100.0 * q, 100.0 * (q + 1.0)))
    assert s3.lag_window(s3.System("y", "f", 30.0, 1.0, 0.5, 2.0))["lead"] == "B"


def test_template_is_not_interpolated_across_gaps():
    t = np.array([0.0, 1.0, 2.0, 100.0, 101.0])
    f = np.arange(5.0)
    out = s3.interp_template(t, f, np.array([0.5, 50.0, 100.0, 200.0]), max_gap=10.0)
    assert out[0] == 0.5 and np.isnan(out[1]) and out[2] == 3.0 and np.isnan(out[3])


def _synthetic(r_true, lag, seed=1):
    rng = np.random.default_rng(seed)
    grid = np.arange(-400.0, 2400.0, 1.0)
    # fast variability (damped random walk, tau = 20 d): a copy at 150 d is not collinear
    s = np.zeros(grid.size)
    for k in range(1, grid.size):
        s[k] = s[k - 1] * np.exp(-1 / 20) + rng.normal(0, 0.05)
    t = np.sort(rng.choice(np.arange(0.0, 2000.0), 500, replace=False))
    tau = 30.0
    fa = 1 + np.interp(t, grid, s)
    fb = 1 + np.interp(t - tau, grid, s) + r_true * np.interp(t - tau - lag, grid, s)
    e = np.full(t.size, 0.002)
    return t, fa + rng.normal(0, 0.002, t.size), fb + rng.normal(0, 0.002, t.size), e, tau


def test_fit_copy_recovers_an_injected_copy_and_gives_zero_without_one():
    p = s3.Params(max_gap=10.0)
    t, fa, fb, e, tau = _synthetic(0.1, 150.0)
    r, m, n = s3.fit_copy(t, fb, e, t, fa, tau, 150.0, p)
    assert n > 200 and abs(m - 1) < 0.1 and abs(r - 0.1) < 0.03
    t, fa, fb, e, tau = _synthetic(0.0, 150.0)
    assert abs(s3.fit_copy(t, fb, e, t, fa, tau, 150.0, p)[0]) < 0.03


def test_window_stat_ignores_nan_and_empty_windows():
    lags = np.array([0.0, 1.0, 2.0])
    assert s3.window_stat(lags, np.array([0.1, np.nan, 0.3]), 0.0, 2.0) == 0.3
    assert np.isnan(s3.window_stat(lags, np.array([0.1, 0.2, 0.3]), 5.0, 6.0))
