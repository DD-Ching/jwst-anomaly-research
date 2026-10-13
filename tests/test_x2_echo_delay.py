"""Offline tests for scripts/x2_echo_delay.py."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest
from scipy.integrate import quad

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import x2_echo_delay as xe  # noqa: E402


@pytest.mark.parametrize("a", [2.05, 2.3, 2.7])
def test_split_integral_matches_direct_quadrature(a):
    rp = math.sqrt(9 - a * a)
    direct = quad(lambda r: 1 / (1 - 2 / math.sqrt(r * r + a * a)), 0, rp, limit=400)[0]
    assert abs(xe.half_cavity(a) - direct) / direct < 1e-8


@pytest.mark.parametrize("d", [1e-8, 1e-6, 1e-4])
def test_small_eps_power_law(d):
    a = 2 + d
    eps = 1 - 2 / a
    approx = 2 * math.pi * math.sqrt(a**3 / eps)
    assert abs(xe.echo_delay(a) - approx) / approx < 0.02


def test_inverse_and_monotonic():
    assert xe.echo_delay(2.001) > xe.echo_delay(2.01) > xe.echo_delay(2.5)
    d = xe.a_for_delay(500.0)
    assert abs(xe.echo_delay(2 + d) - 500.0) < 1e-6


def test_requires_two_peaks():
    with pytest.raises(ValueError):
        xe.half_cavity(3.2)
