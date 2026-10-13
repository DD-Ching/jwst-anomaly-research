"""Fast offline tests for scripts/x2_scatter.py (X2 S1 tool)."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import x2_scatter as xs  # noqa: E402


def test_em_peak_is_photon_sphere_scalar_is_not():
    _, R1, v1 = xs.peak(2.001, 2, 1)
    _, R0, _ = xs.peak(2.001, 2, 0)
    assert abs(R1 - 3.0) < 1e-6
    assert abs(v1 - 6.0 / 27.0) < 1e-9
    assert abs(R0 - 3.0) > 0.01  # the scalar peak is not at the photon sphere


def test_rw_proxy_barrier_height_equals_schwarzschild_spin2():
    # the proxy is the Schwarzschild l = 2 Regge-Wheeler form in R, so its peak is the same 0.151287 / M^2
    _, _, v = xs.peak(2.0005, 2, "rw")
    assert abs(v - 0.151287) < 1e-5


@pytest.mark.parametrize("s", [0, 1])
def test_flux_and_two_routes(s):
    g = xs.Grid(2.01, 2, s, xmax=400.0, h=0.01)
    _, _, vp = xs.peak(2.01, 2, s)
    w = math.sqrt(vp)
    R2, T2 = xs.single_barrier(g, np.array([w]))
    assert abs(R2[0] + T2[0] - 1) < 1e-8
    _, t2b = xs.single_barrier_ode(g, w)
    assert abs(T2[0] - t2b) / T2[0] < 1e-4


def test_schwarzschild_low_frequency_absorption():
    w = 0.002
    ratio = xs.schwarzschild_scalar_T2(np.array([w]))[0] / (16 * w * w)
    assert abs(ratio - 1) < 0.02  # -> 1 as w -> 0 (1.013 at w = 0.002)
