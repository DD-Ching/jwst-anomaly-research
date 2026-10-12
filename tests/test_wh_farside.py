"""Offline tests for scripts/wh_farside.py."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import wh_farside as wf  # noqa: E402


@pytest.mark.parametrize("a", [2.1, 3.0, 4.0, 6.0, 9.0])
def test_b_crit_matches_numerical_barrier(a):
    assert abs(wf.b_crit(a) - wf.b_crit_numeric(a)) < 1e-3


def test_theta_g_matches_gravity_value():
    # GRAVITY 2023: R_g = 5.12 uas for M = 4.297e6 Msun at R0 = 8.277 kpc
    assert abs(wf.theta_g(4.297e6) - 5.12) < 0.01


def test_kepler_radius_and_period_roundtrip():
    # P = 2 pi R^(3/2) GM/c^3 for Omega^2 = M/R^3
    m = 4.297e6
    r = wf.kepler_radius(60.0, m)
    tg = wf.G * m * wf.MSUN / wf.C**3
    assert abs(2 * math.pi * r**1.5 * tg / 60.0 - 60.0) < 1e-9
    assert 8.9 < r < 9.1


def test_eht_shadow_delta_zero_below_photon_sphere():
    assert wf.shadow_delta(2.5) == 0.0
    assert wf.shadow_delta(5.0) > 0.2


def test_mass_estimator_returns_input_mass_for_same_side_loop():
    # same-side loop r_o = (R_e + 1) theta_g must give back M (Appendix A of GRAVITY 2023)
    m = 4.297e6
    ro = (wf.kepler_radius(60.0, m) + 1.0) * wf.theta_g(m)
    assert abs(wf.gravity_mass_estimator(ro, 60.0) / m - 1) < 1e-3
