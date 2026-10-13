"""Offline tests for scripts/x1_s3.py (X1 S3 fit), on synthetic data only."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import x1_s2  # noqa: E402
import x1_s3  # noqa: E402

FAR_WARN = pytest.mark.filterwarnings("ignore::scipy.integrate.IntegrationWarning")


def _synthetic(track, period, theta, phi0, sense, centre):
    t = np.linspace(0.0, 30.0, 10)
    cs = track[:, 1:3] * x1_s2.THETA_G
    arg = (phi0 + 2 * math.pi * t / period) % (2 * math.pi)
    pp = np.concatenate([track[:, 0], track[:, 0] + 2 * math.pi])
    spot = np.stack([np.interp(arg, pp, np.tile(cs[:, k], 2)) for k in (0, 1)], -1)
    spot = spot * np.array([1.0, sense])
    rot = np.array([[math.cos(theta), -math.sin(theta)], [math.sin(theta), math.cos(theta)]])
    xy = spot @ rot.T + np.asarray(centre)
    return {"t": t, "xy": xy, "err": np.full_like(xy, 10.0), "flux": np.linspace(0.3, 0.6, 10)}


@FAR_WARN
@pytest.mark.parametrize("a,side", [(0.0, "near"), (5.0, "far")])
def test_fit_recovers_injected_track(a, side):
    tr, _ = x1_s2.track(a, 9.0, side, 1.0)
    g = 2 * math.pi / x1_s3.N_GRID
    d = _synthetic(tr, 60.0, 10 * g, 30 * g, -1.0, (40.0, -25.0))
    best = x1_s3.best_chi2(tr, d, 60.0)
    assert best["chi2"] < 1e-6
    assert best["sense"] == -1.0
    assert best["centre_uas"] == pytest.approx([40.0, -25.0], abs=1e-6)


@FAR_WARN
def test_far_side_cannot_fit_an_ordinary_loop():
    near, _ = x1_s2.track(0.0, 9.0, "near", 1.0)
    far, _ = x1_s2.track(3.0, 9.0, "far", 1.0)
    d = _synthetic(near, 60.0, 0.0, 0.0, 1.0, (0.0, 0.0))
    assert x1_s3.best_chi2(far, d, 60.0)["chi2"] > 25 * x1_s3.best_chi2(near, d, 60.0)["chi2"] + 50
