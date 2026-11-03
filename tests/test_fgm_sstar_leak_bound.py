"""Offline tests for scripts/fgm_sstar_leak_bound.py (cheap checks plus the stored bound)."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import fgm_sstar_leak_bound as sb  # noqa: E402

STORED = json.loads((sb.OUT / "bound.json").read_text())


def test_thresholds_reproduced():
    t = STORED["thresholds"]
    assert math.isclose(t["N_th"], 0.2554, abs_tol=5e-5)
    assert math.isclose(t["closure"], 0.3125)
    assert math.isclose(t["K_critical"], 1.2235, abs_tol=5e-4)


def test_rn_geodesics():
    # photon sphere and the weak-field deflection of RN (second order)
    assert math.isclose(sb.f(sb.R_PH) * sb.B_C**2, sb.R_PH**2, rel_tol=1e-12)
    b = 50.0
    defl = 2 * sb.leg_sweep(1e12, b) - math.pi
    second = 4 * sb.M / b + (15 * math.pi / 4) * sb.M**2 / b**2 - (3 * math.pi / 4) * sb.Q2 / b**2
    assert abs(defl - second) < 3e-4


def test_edge_saddle_formula():
    for q in (100.0, 300.0):
        assert abs(sb.envelope_ratio(q, span=24) - 1) < 0.02


def test_stored_bound_is_consistent():
    g = STORED["geometry"]
    assert g["first_orbit_after_family2"] > 15
    assert STORED["aperture_scan"]["min_family2_sigma_over_omega"] > 0.4
    assert STORED["aperture_scan"]["earliest_orbit_minus_family2_edge"] > 15
    e = STORED["edge"]
    assert e["family1_over_family2_amplitude"] < 1e-4
    assert abs(e["envelope_ratio_mp_at_q2"] - 1) < 0.02
    b = STORED["bound"]
    assert b["log10_K_lower"] > 30 and STORED["verdict"].startswith("exclusion supported")
