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


def test_null_energy_is_violated_only_in_the_shortcut_direction():
    for eta in (1e-2, 0.95):
        nc = ks.null_and_causal(eta)
        assert nc["outgoing_Rll_symbolic"] == "0"  # exterior light direction: NEC saturated
        assert nc["outgoing_fixed_rho_radial_acceleration"] == "0"  # those curves are geodesics
        assert nc["incoming_Rll_min_eps1"] < 0  # shortcut direction: NEC violated in the wall
        assert nc["t_is_time_function"] and nc["k_min"] == pytest.approx(1 - eta, rel=1e-6)


def test_curved_qi_corrections_small_for_latency_rows_only():
    near_one = ks.kontou_olum_terms(1 - 1e-9)
    assert near_one["t0_over_eps"] < 0.01  # sampling far inside the curvature-variation scale
    for key in ("ricci_term_rel", "second_derivative_term_rel", "third_derivative_term_rel"):
        assert near_one[key] < 1e-3
    one_second_advance = ks.kontou_olum_terms(1.6e-9)
    assert one_second_advance["t0_over_eps"] > 10  # tau0 = sigma r_c spans many wall thicknesses
    assert one_second_advance["third_derivative_term_rel"] > 1  # first-order expansion fails
