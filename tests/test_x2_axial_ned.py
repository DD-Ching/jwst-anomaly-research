"""Fast offline tests for scripts/x2_axial_ned.py (X2 S4)."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import x2_axial_ned as an  # noqa: E402


def _d(fun, r, h=1e-4):
    return (fun(r + h) - fun(r - h)) / (2 * h)


def test_v_af_equals_its_definition_and_regge_wheeler_limit():
    a, r = 2.3, np.array([0.0, 0.7, 1.9, 4.2, 9.0])
    R = lambda q: np.sqrt(q * q + a * a)  # noqa: E731
    f = lambda q: 1 - 2 / R(q)  # noqa: E731
    # R (1/R)_{,r*r*} = R f d/dr ( f d/dr (1/R) )
    inner = lambda q: f(q) * _d(lambda s: 1 / R(s), q)  # noqa: E731
    want = f(r) * an.MU2 / R(r) ** 2 + R(r) * f(r) * _d(inner, r)
    assert np.allclose(an.v_af(r, a), want, rtol=1e-6, atol=1e-10)
    rr = np.array([3.0, 5.0])
    rw = (1 - 2 / rr) * (6 / rr**2 - 6 / rr**3)
    assert np.allclose(an.v_af(rr, 1e-7), rw, rtol=1e-9)


def test_v11_v12_equal_definitions_through_LF():
    # L_F proportional to 1 / (R^4 (G00 - G22)), G00 - G22 = 3 a^2 / R^5  ->  L_F ~ R
    a, q = 2.4, 0.8
    r = np.array([0.0, 0.5, 1.7, 4.0])
    R = lambda s: np.sqrt(s * s + a * a)  # noqa: E731
    f = lambda s: 1 - 2 / R(s)  # noqa: E731
    LF = lambda s: 8 * math.pi * q * q / (R(s) ** 4 * 3 * a * a / R(s) ** 5)  # noqa: E731
    inner = lambda s: f(s) / np.sqrt(LF(s)) * _d(LF, s)  # noqa: E731
    v11 = (
        f(r)
        / R(r) ** 2
        * (
            an.MU2
            + 2
            + R(r) ** 2 / (2 * np.sqrt(LF(r))) * _d(inner, r)
            + 16 * math.pi * q * q / (R(r) ** 2 * LF(r))
        )
    )
    v12 = 2 * math.sqrt(an.MU2) * q * f(r) / (R(r) ** 3 * np.sqrt(LF(r)))
    assert np.allclose(an.v11(r, a), v11, rtol=1e-6)
    assert np.allclose(an.v12(r, a), v12, rtol=1e-9)


def test_coupling_off_gw_channel_equals_single_channel():
    a = 2.2

    def init(x, dt):
        g = np.exp(-((x - 5.0) ** 2) / 4.0)
        return [(g, g.copy())]

    _, r_off = an.evolve("ned_off", a, 150.0, 150.0, 0.2, 100.0, init, [30.0])
    _, r_af = an.evolve("af", a, 150.0, 150.0, 0.2, 100.0, init, [30.0])
    assert np.allclose(r_off[0, 0], r_af[0, 0], atol=1e-14)
    assert np.max(np.abs(r_off[0, 1])) == 0.0


def test_matrix_pencil_recovers_damped_modes():
    t = np.arange(0, 400, 1.0)
    y = np.exp(-0.01 * t) * np.cos(0.3 * t) + 0.3 * np.exp(-0.002 * t) * np.cos(0.21 * t + 1)
    w, amp = an.matrix_pencil(y, 1.0, order=4)
    w = w[np.argsort(-amp)]
    assert abs(w[0] - complex(0.3, -0.01)) < 1e-6
    assert abs(w[1] - complex(0.21, -0.002)) < 1e-6


def test_gw_peak_far_from_throat_near_a_2():
    # a - 2 = 3.6e-4: the tortoise distance from the throat to the peak is ~ 330 M (dt/4)
    x = an.gw_peak_x("af", 2.000358)
    assert 250.0 < x < 450.0
