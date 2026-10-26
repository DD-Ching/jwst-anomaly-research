"""Offline tests of scripts/w3_shortcad_reach.py (results/w3_kepler/shortcad/PLAN.md)."""

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import w3_shortcad_reach as ws  # noqa: E402


def test_sc_months_reads_one_digit_per_quarter():
    assert ws.sc_months("kplr008264588_sc_Q010000000000000000") == (1, ws.SC_MONTH_D)
    n, d = ws.sc_months("kplr000000001_sc_Q113000000000000002")
    assert n == 7 and d == ws.SC_Q0_D + 6 * ws.SC_MONTH_D


def test_identifiable_needs_dark_run_and_both_shoulders():
    dark = [0.0, 0.0, 0.0]
    assert ws.identifiable(np.array([1, 1.2, 1.2] + dark + [1.2, 1.2, 1])) == (True, True)
    assert ws.identifiable(np.array([1, 1.2, 1.0] + dark + [1.2, 1.2, 1])) == (True, False)
    assert ws.identifiable(np.array([1.2, 1.2, 0, 0, 1, 0, 1.2, 1.2])) == (False, False)


def test_level0_matches_the_analytic_rate():
    # w = 2: N = T v (pi / 2) sqrt(4 G M_sun / c^2) D^1.5 M^-0.5 (int sqrt(x(1-x)) dx = pi / 8)
    d, t, m, v = 500.0, 2.0, 1e-3, 220.0
    n = ws.events_per_density([d], [1.0], [t], m, v, 60.0)
    v_pc_yr = v * 365.25 * 86400 / ws.PC_KM
    expect = t * v_pc_yr * math.pi / 2 * math.sqrt(ws.SCHW4_PC) * d**1.5 / math.sqrt(m)
    assert abs(n / expect - 1) < 2e-3


def test_level1_never_exceeds_level0_and_vanishes_for_huge_rho():
    lr, lt = np.log10(ws.RHO_GRID), np.log10(ws.TAU_GRID)
    shape = (lr, lt, np.full((len(lr), len(lt)), 1.9))
    n0 = ws.events_per_density([300.0], [1.0], [1.0], 1e-6, 30.0, 20.0)
    n1 = ws.events_per_density([300.0], [1.0], [1.0], 1e-6, 30.0, 20.0, shape)
    assert 0 <= n1 <= n0
    assert ws.w_lookup(np.array([10.0]), np.array([100.0]), shape)[0] == 0.0
    assert ws.w_lookup(np.array([0.1]), np.array([1e6]), shape)[0] == 1.9  # clamped above the grid


def test_gate():
    assert ws.gate([1, 12, 15, None, 3])["passes"]
    assert ws.gate([1, 12, 15, None, 3])["intervals_idx"] == [[1, 3]]
    assert not ws.gate([12, 15, 2, 30, 40])["passes"]
