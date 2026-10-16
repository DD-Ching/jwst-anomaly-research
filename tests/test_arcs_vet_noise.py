"""Tests for scripts/arcs_vet_noise.py (synthetic values; no network)."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest

_DIR = Path(__file__).resolve().parents[1] / "scripts"
_spec = importlib.util.spec_from_file_location("arcs_vet_noise", _DIR / "arcs_vet_noise.py")
vn = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(vn)


def test_gate_ok():
    assert vn.gate_ok([0.95, 1.0])
    assert not vn.gate_ok([0.94, 1.0])
    assert not vn.gate_ok([])
    assert not vn.gate_ok([np.nan, 1.0])


def test_noise_class():
    assert vn.noise_class(0.05) == "noise-explicable"
    assert vn.noise_class(0.049) == "not noise-explicable"
    assert vn.noise_class(float("nan")) == "undetermined"


def test_excess_p_matches_binomial():
    assert vn.excess_p([0.5, 0.5], 1) == pytest.approx(0.75)
    assert vn.excess_p([0.5, 0.5], 0) == pytest.approx(1.0)
    assert vn.excess_p([0.5, 0.5], 3) == 0.0
    assert vn.excess_class(0.05) == "excess explained by measurement noise"
    assert vn.excess_class(0.01) == "excess remains"


@pytest.mark.parametrize("theta_deg", [0.0, 35.0, 120.0])
def test_noiseless_injection_recovers_axis(theta_deg):
    shape = (81, 81)
    stamp = vn.gaussian_stamp(shape, 40, 40, 100.0, 3.0, 0.5, np.radians(theta_deg))
    assert stamp.sum() == pytest.approx(100.0, rel=1e-3)
    t, e = vn.measure(stamp, 33, 3, 3.0)
    assert vn.abs_.axis_offset(t, theta_deg) < 2.0
    assert 0.0 < e < 0.5
