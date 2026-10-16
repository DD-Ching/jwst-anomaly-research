"""Tests for scripts/arcs_faint.py (synthetic values; no network)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("arcs_faint", _DIR / "arcs_faint.py")
af = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(af)


def test_faint_rank_p_all_anti_faintest():
    mag = [20.0, 21.0, 22.0, 23.0, 28.0, 27.0]
    anti = [False, False, False, False, True, True]
    r = af.faint_rank_p(mag, anti)
    assert r["k_faintest"] == 2
    assert r["p"] == 1 / 15  # C(2,2) / C(6,2)


def test_faint_rank_p_brightest_anti_sets_k():
    mag = [20.0, 28.0, 27.0, 26.0]
    anti = [True, True, False, False]
    r = af.faint_rank_p(mag, anti)
    assert r["k_faintest"] == 4 and r["p"] == 1.0


def test_faint_rank_p_no_anti():
    assert af.faint_rank_p([20.0, 21.0], [False, False])["p"] == 1.0


_s2 = importlib.util.spec_from_file_location("arcs_second_models", _DIR / "arcs_second_models.py")
asm = importlib.util.module_from_spec(_s2)
_s2.loader.exec_module(asm)
_s3 = importlib.util.spec_from_file_location("arcs_bgsub_shapes", _DIR / "arcs_bgsub_shapes.py")
abs_ = importlib.util.module_from_spec(_s3)
_s3.loader.exec_module(abs_)


def test_stretch_angle_flips_at_kappa_one():
    import numpy as np

    a = asm.stretch_angle_px(np.array([0.5, 1.5]), np.array([0.2, 0.2]), np.array([0.0, 0.0]))
    assert a[0] == 0.0 and a[1] == 90.0
    assert asm.axis_offset(np.array([170.0]), np.array([10.0]))[0] == 20.0


def test_weighted_axis_on_gradient():
    import numpy as np

    yy, xx = np.indices((41, 41), dtype=float)
    dx, dy = xx - 20, yy - 20
    th = np.radians(30.0)
    u, v = dx * np.cos(th) + dy * np.sin(th), -dx * np.sin(th) + dy * np.cos(th)
    img = np.exp(-0.5 * (u**2 / 9.0 + v**2 / 1.0))
    theta, e = abs_.weighted_axis(img, 20.0, 20.0, 3.0)
    assert abs(theta - 30.0) < 1.0 and 0.2 < e < 0.7
