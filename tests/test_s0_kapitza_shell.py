"""Tests for scripts/s0_kapitza_shell.py (thin-shell throat with a confined radiation component)."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np

_P = Path(__file__).resolve().parents[1] / "scripts" / "s0_kapitza_shell.py"
_spec = importlib.util.spec_from_file_location("s0_kapitza_shell", _P)
k = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(k)


def test_static_throat_is_an_equilibrium():
    for a0 in (2.5, 4.0, 10.0):
        v0, v1, _ = k.derivatives(a0, 0.3, 0.5)
        assert abs(v0) < 1e-12 and abs(v1) < 1e-5


def test_poisson_visser_beta_slope():
    # without radiation V''(a0) is linear in beta^2 with slope -(4/a0^2)(1 - 3M/a0)
    for a0 in (4.0, 6.0):
        s = (k.derivatives(a0, 1.0, 0.0)[2] - k.derivatives(a0, 0.0, 0.0)[2]) / 1.0
        assert abs(s + 4 / a0**2 * (1 - 3 / a0)) < 1e-5


def test_radiation_effect_sign_follows_effective_eos():
    # for a0 > 3M the shell has sigma + p < 0, and radiation pushes beta2_eff away from 1/2
    a0 = 6.0
    assert np.sign(k.static(a0)[0] + k.static(a0)[1]) < 0
    assert k.derivatives(a0, 0.0, 0.1)[2] > k.derivatives(a0, 0.0, 0.0)[2]
    assert k.derivatives(a0, 1.0, 0.1)[2] < k.derivatives(a0, 1.0, 0.0)[2]
