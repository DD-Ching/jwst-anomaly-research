"""Offline tests for scripts/goal_krasnikov_source.py (Krasnikov tube source checks)."""

import importlib.util
import sys
from pathlib import Path

import pytest

_DIR = Path(__file__).resolve().parents[1] / "scripts"
_spec = importlib.util.spec_from_file_location(
    "goal_krasnikov_source", _DIR / "goal_krasnikov_source.py"
)
ks = importlib.util.module_from_spec(_spec)
sys.modules["goal_krasnikov_source"] = ks
_spec.loader.exec_module(ks)


def test_weyl_is_traceless_and_below_ricci_on_the_wall():
    for eta in (1e-6, 1.0):
        w = ks.weyl_ricci(eta)
        assert w["weyl_trace_residual"] < 1e-9 * max(1.0, w["max_abs_ricci"])
        assert w["ratio_max"] < 0.5  # Krasnikov 2003 Sec. III A loophole needs ratio >> 1


def test_riemann_reproduces_round5_curvature_coefficient():
    import numpy as np

    r, k, k1, k2 = ks.fc.wall_profile(1e-6, 1.0)
    rmax = np.max(np.abs(ks.riemann_array(r, k, k1, k2)))
    assert rmax / 1e-6 == pytest.approx(3.0792, rel=1e-3)


def test_latency_target_needs_no_sampling_extrapolation():
    rows = ks.latency_rows()
    one_second = [r for r in rows if r["latency_s"] == 1.0]
    for r in one_second:
        assert r["rc_over_eps"] < 1  # tau0 = sigma r_c is below sigma eps: inside ER's own rule
        assert r["eps_max_planck_tau_sigma_rc"] == r["eps_max_planck_tau_sigma_min_rc_eps"]
        assert 5e3 < r["eps_max_planck_tau_sigma_rc"] < 2e4


def test_integrated_negative_energy_and_self_consistency():
    w = ks.wall_energy_integrals(1.0)
    assert w["I_neg"] < 0 < w["I_pos"] and 0.2 < w["effective_alpha"] < 0.5
    sc = ks.self_consistency(1 - 1e-9)
    assert sc["eps_self_consistent_planck_N1"] < 1.0  # one field's polarisation: sub-Planck wall
    assert sc["fields_needed_to_reach_qi_thickness"] == pytest.approx(270 / ks.P.sigma**4, rel=1e-6)
