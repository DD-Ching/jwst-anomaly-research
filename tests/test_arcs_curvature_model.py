"""Tests for scripts/arcs_curvature_model.py (an SIS stand-in model; no network)."""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import numpy as np

_DIR = Path(__file__).resolve().parents[1] / "scripts"
_spec = importlib.util.spec_from_file_location(
    "arcs_curvature_model", _DIR / "arcs_curvature_model.py"
)
acm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(acm)


class _SIS:
    """Singular isothermal sphere at the frame origin; sky = frame (x West, y North, arcsec)."""

    theta_e = 10.0

    def to_frame(self, ra, dec):
        return np.asarray(ra, float), np.asarray(dec, float)

    def source_points(self, x, y, z_s):
        r = np.hypot(x, y)
        return x - self.theta_e * x / r, y - self.theta_e * y / r


def test_predicted_arc_bends_around_the_lens():
    args = argparse.Namespace(sigma_src=0.1, size=8.0, step=0.03, frac=0.1)
    # image just outside the Einstein ring, east of the lens (x = -10.5 is East)
    fit = acm.predicted_fit(_SIS(), -10.5, 0.0, 2.0, args)
    # the lens is to the west (+x West, so east component -1): curvature side points west
    assert fit["side_e"] < -0.9
    assert abs(abs(fit["c"]) - 0.5 / 10.5) < 0.02
    assert fit["length"] > 2.0
