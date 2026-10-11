"""Offline tests for scripts/goal_ea_response.py (Einstein-aether spin-0 channel)."""

import importlib.util
import sys
from pathlib import Path

import sympy as sp

_DIR = Path(__file__).resolve().parents[1] / "scripts"
_spec = importlib.util.spec_from_file_location("goal_ea_response", _DIR / "goal_ea_response.py")
ea = importlib.util.module_from_spec(_spec)
sys.modules["goal_ea_response"] = ea
_spec.loader.exec_module(ea)


def _zero(expr: str) -> bool:
    return sp.simplify(sp.sympify(expr)) == 0


def _equal(expr: str, expected: str) -> bool:
    return sp.simplify(sp.sympify(expr) - sp.sympify(expected)) == 0


def test_positive_controls_speeds_newton_and_gr_quadrupole():
    c = ea.controls()
    for mode in ("S", "V", "T"):
        assert len(c["pole_minus_literature_speed2"][mode]) == 1
        assert _zero(c["pole_minus_literature_speed2"][mode][0])
    assert _zero(c["static_Ezz_minus_4piGN_rho"])
    assert _equal(c["tensor_kernel_over_GR"], "1")


def test_spherical_source_exterior_is_causal_spin0_only():
    s = ea.spherical_source()
    assert s["no_instantaneous_term"]
    assert _equal(s["pole_Exx_over_alpha_cS2_k4_Q"], "1")
    assert _equal(s["kappa_leading_c14_to_0"], "c14/4")


def test_far_field_breathing_matches_yagi_eq106():
    ff = ea.far_field_general()
    assert _equal(ff["xx"]["leading_c14_to_0"], "-c14*(Ixx + Iyy + Izz)/4")
    y = ea.yagi_comparison()
    assert _zero(y["Ikk_ours_minus_yagi"])
    assert _equal(y["Qnn_ours_over_yagi"], "1")


def test_general_source_kernel_matches_all_components():
    for v in ea.general_source_kernel().values():
        assert _zero(v)


def test_vector_sector_bound():
    v = ea.vector_sector_c13()
    assert v["max_kappaV_times_cV2_(1-c13)_on_grid"] <= 1.0 + 1e-9


def test_link_budget_is_hopeless_and_worsens_with_speed():
    rows = ea.link_budget()
    for r in rows:
        assert r["integration_time_s_for_snr"] > 1e40  # vs 4e17 s, the age of the universe
    one_s = [r for r in rows if r["D_ly"] == 20.0 and r["latency_s"] == 1.0][0]
    one_yr = [r for r in rows if r["D_ly"] == 20.0 and r["latency_s"] > 1e7][0]
    assert one_s["kappa"] < 1e-18 < one_yr["kappa"]
