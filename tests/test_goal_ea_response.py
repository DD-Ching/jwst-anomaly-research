"""Offline tests for scripts/goal_ea_response.py (Einstein-aether spin-0 channel)."""

import importlib.util
import sys
from pathlib import Path

_DIR = Path(__file__).resolve().parents[1] / "scripts"
_spec = importlib.util.spec_from_file_location("goal_ea_response", _DIR / "goal_ea_response.py")
ea = importlib.util.module_from_spec(_spec)
sys.modules["goal_ea_response"] = ea
_spec.loader.exec_module(ea)


def test_positive_controls_speeds_newton_and_gr_quadrupole():
    c = ea.controls()
    for mode in ("S", "V", "T"):
        assert c["pole_minus_literature_speed2"][mode] == ["0"]
    assert c["static_Ezz_minus_4piGN_rho"] == "0"
    assert c["tensor_kernel_over_GR"] == "1"


def test_spherical_source_exterior_is_causal_spin0_only():
    s = ea.spherical_source()
    assert s["no_instantaneous_term"]
    assert s["kappa_leading_c14_to_0"] == "c14/4"


def test_far_field_breathing_scales_with_c14_and_matches_yagi_on_ikk():
    ff = ea.far_field_general()
    assert ff["xx"]["leading_c14_to_0"] == "-c14*(Ixx + Iyy + Izz)/4"
    y = ea.yagi_comparison()
    assert y["Ikk_ours_minus_yagi"] == "0"
    assert y["Qnn_ours_over_yagi"] == "1/2"  # unresolved factor at O(c14^2), recorded


def test_link_budget_is_hopeless_and_worsens_with_speed():
    rows = ea.link_budget()
    for r in rows:
        assert r["integration_time_s_for_snr"] > 1e40  # vs 4e17 s, the age of the universe
    one_s = [r for r in rows if r["D_ly"] == 20.0 and r["latency_s"] == 1.0][0]
    one_yr = [r for r in rows if r["D_ly"] == 20.0 and r["latency_s"] > 1e7][0]
    assert one_s["kappa"] < 1e-18 < one_yr["kappa"]
