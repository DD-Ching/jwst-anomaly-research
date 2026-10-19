"""Offline tests for scripts/smacs_at2024jfz.py: the frozen rule and the delay integral (SIS)."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import smacs_at2024jfz as s  # noqa: E402


def test_testable_window_and_flux():
    t0 = 100.0
    epochs = (90.0, 120.0, 400.0)
    # arrives at 110: epoch 120 is +10 d after arrival (inside), 90 is -20 d (outside)
    assert s.testable(10.0, 1.0, t0=t0, epochs=epochs) == [120.0]
    # arrives at 122: epoch 120 is -2 d (edge, inside)
    assert s.testable(22.0, 1.0, t0=t0, epochs=epochs) == [120.0]
    assert s.testable(22.1, 1.0, t0=t0, epochs=epochs) == []
    # arrives at 370: epoch 400 is +30 d (edge, inside)
    assert s.testable(270.0, 1.0, t0=t0, epochs=epochs) == [400.0]
    # too faint: 2579 * 0.05 < 170 nJy
    assert s.testable(10.0, 0.05, t0=t0, epochs=epochs) == []
    assert s.testable(10.0, -0.07, t0=t0, epochs=epochs) == [120.0]  # |mu ratio| counts


def test_sis_delay_matches_analytic():
    # SIS with theta_E = 1": alpha = theta_E * theta / |theta|, psi = theta_E |theta|.
    te = 1.0

    def deflect(x, y):
        r = np.hypot(x, y)
        return te * x / r, te * y / r

    beta = 0.3
    p_plus, p_minus = (beta + te, 0.0), (beta - te, 0.0)
    # path that avoids the singular centre: via a corner off-axis
    dpsi = s.psi_difference(deflect, p_plus, p_minus, via=(0.0, 1.0))
    # straight path through the corner pieces is not used here; compare with the analytic psi
    assert dpsi == pytest.approx(te * abs(p_minus[0]) - te * abs(p_plus[0]), abs=2e-3)
    phi = [0.5 * (p[0] - beta) ** 2 for p in (p_plus, p_minus)]
    dphi = (phi[1] - dpsi) - phi[0]
    assert dphi == pytest.approx(
        2 * te * beta, abs=2e-3
    )  # analytic SIS: delta phi = 2 theta_E beta


def test_path_independence_smooth_field():
    # alpha = grad psi for psi = x^2 y + sin(y): any path gives the same difference
    def deflect(x, y):
        return 2 * x * y, x**2 + np.cos(y)

    p0, p1 = (0.2, -0.5), (1.3, 0.8)
    exact = (p1[0] ** 2 * p1[1] + np.sin(p1[1])) - (p0[0] ** 2 * p0[1] + np.sin(p0[1]))
    assert s.psi_difference(deflect, p0, p1) == pytest.approx(exact, rel=1e-9)
    assert s.psi_difference(deflect, p0, p1, via=(p1[0], p0[1])) == pytest.approx(exact, rel=1e-9)
