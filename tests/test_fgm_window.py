"""Offline tests for scripts/fgm_window.py (analytic FGM window, a few seconds)."""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import fgm_window as fw  # noqa: E402


def test_green_function_closed_form_matches_series():
    # FGM (3.14) against the harmonic sum (3.15)
    for th, tol in ((math.pi / 2, 1e-5), (math.pi, 1e-3)):
        assert math.isclose(fw.green_h(th), fw.green_h_series(th), rel_tol=tol)


def test_window_reproduces_fgm_transit_time():
    # t_min from DeltaV equals FGM (3.17) read with H / r+^2, for two charge ratios
    for r_minus in (0.5, 0.2):
        w = fw.window(math.pi / 2, 1e-70, 20.0, x0=5e3, r_minus=r_minus)
        t317 = fw.t_min_fgm_317(math.pi / 2, 1e-70, 20.0, x0=5e3, r_minus=r_minus)
        assert math.isclose(w["t_min"], t317, rel_tol=1e-12)
        # the entry deadline is minus half the minimal transit plus ln 2 / kappa
        k = fw.kappa(r_minus)
        assert math.isclose(w["v_max"], -w["t_min"] / 2 + math.log(2) / k, rel_tol=1e-12)


def test_closure_threshold_blueshifts_into_the_past():
    a, b = fw.de_crit(-100.0, 20.0), fw.de_crit(-104.0, 20.0)
    assert math.isclose(a / b, math.exp(4 * fw.kappa()), rel_tol=1e-12)


def test_packets():
    g = fw.packet_gaussian(500.0, 100.0)
    assert math.isclose(g["sd_v"], 0.005)
    assert g["weight_above_5"] > 0.9999999 and 0 < g["weight_below_1"] < 1e-7
    b = fw.packet_bump(500.0, 200.0, n=4001)
    assert abs(b["mean_E"] - 500.0) < 1e-6 and b["weight_below_1"] == 0.0
    assert b["prob_within_0.5"] > 0.9999
