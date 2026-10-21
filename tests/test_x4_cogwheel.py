"""Offline tests for scripts/x4_cogwheel.py (X4 cogwheel round): the map from cogwheel parameters to
our box, the in-box flag and the periodic box distance. No lalsuite or cogwheel needed."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

_S = Path(__file__).resolve().parents[1] / "scripts"


def _load(name):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _S / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


X = _load("x4_cogwheel")
e = X.e


def cog_par(x):
    """Our box point -> cogwheel parameters, rotating in-plane spins by +phi_ref by hand."""
    mc, q, s1, s2, inc, phi = e.to_cart(np.asarray(x, float))
    m1, m2 = X.b.masses(mc, q)
    c, s = np.cos(phi), np.sin(phi)
    par = {"m1": m1, "m2": m2, "iota": inc, "phi_ref": phi}
    for i, sp in ((1, s1), (2, s2)):
        par[f"s{i}x_n"] = c * sp[0] - s * sp[1]
        par[f"s{i}y_n"] = s * sp[0] + c * sp[1]
        par[f"s{i}z"] = sp[2]
    return par


def test_round_trip_and_inside():
    rng = np.random.default_rng(1)
    lo = np.array([bd[0] for bd in e.BOUNDS])
    hi = np.array([bd[1] for bd in e.BOUNDS])
    for _ in range(20):
        x = lo + rng.random(lo.size) * (hi - lo)
        x[2], x[5] = min(x[2], 0.98), min(x[5], 0.98)
        x[1] = min(x[1], 14.9)
        y, inside, mags = X.to_box(cog_par(x))
        assert X.box_diff(x, y) < 1e-9
        assert inside
        assert np.allclose(mags, [x[2], x[5]])


def test_outside_flags():
    x = np.array([np.log(e.MC0), 3.0, 0.5, 0.2, 1.0, 0.5, -0.3, 2.0, 1.0, 0.5])
    par = cog_par(x)
    scale = 0.995 / 0.5
    for key in ("s1x_n", "s1y_n", "s1z"):
        par[key] *= scale  # |chi1| = 0.995 > 0.99
    y, inside, mags = X.to_box(par)
    assert not inside and mags[0] > 0.99 and y[2] == 0.99
    par = cog_par(x)
    m_tot = par["m1"] + par["m2"]
    par["m1"], par["m2"] = m_tot * 16 / 17, m_tot / 17  # q = 16 > 15
    assert not X.to_box(par)[1]


def test_box_diff_wraps_angles():
    a = np.zeros(10)
    b = np.zeros(10)
    a[4], b[4] = 0.01, 2 * np.pi - 0.01
    assert abs(X.box_diff(a, b) - 0.02) < 1e-12
    b[0] = 0.5
    assert X.box_diff(a, b) == 0.5
