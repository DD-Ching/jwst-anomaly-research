"""Offline tests for scripts/goal_ftl_candidates.py (Krasnikov tube QI trade-off, Scharnhorst)."""

import importlib.util
import math
import sys
from pathlib import Path

import pytest

_DIR = Path(__file__).resolve().parents[1] / "scripts"
_spec = importlib.util.spec_from_file_location(
    "goal_ftl_candidates", _DIR / "goal_ftl_candidates.py"
)
fc = importlib.util.module_from_spec(_spec)
sys.modules["goal_ftl_candidates"] = fc
_spec.loader.exec_module(fc)


def test_positive_controls_reproduce_everett_roman():
    c = fc.krasnikov_controls()
    assert c["max_abs_diff_vs_ER_eq14_delta_small"] < 1e-8  # independent Einstein tensor = Eq. 14
    assert c["endcap_T_at_rho_max_minus_eps"] == pytest.approx(c["ER_eq38_value"], rel=5e-3)


def test_small_eta_scaling_is_linear():
    c = fc.krasnikov_controls()
    t, r = c["small_eta_T_coeff_over_eta"], c["small_eta_R_coeff_over_eta"]
    assert t[-1] == pytest.approx(t[-2], rel=1e-3) and r[-1] == pytest.approx(r[-2], rel=1e-3)


def test_qi_walls_are_planck_scale_for_the_20_ly_task():
    rows = [r for r in fc.krasnikov_task()["rows"] if r["D_ly"] == 20.0]
    assert all(r["eps_max_planck"] < 2e4 for r in rows)
    one_second = next(r for r in rows if r["advance_s"] == 1.0)
    assert one_second["eps_max_planck"] < 1.0  # sub-Planck: no semiclassical tube at all
    assert all(r["eps_max_loose_m"] < 1e-26 for r in rows)  # looser prescription still microscopic


def test_scharnhorst_advance_is_unresolvable_in_its_own_regime():
    s = fc.scharnhorst()
    assert all(r["advance_times_mec2_over_hbar"] < 1e-6 for r in s["rows"])
    xi = s["xi_1um"]
    want = 11 * math.pi**2 * fc.ALPHA**2 / 4050 / (1e-6 / fc.LAMBDA_C_BAR_M) ** 4
    assert xi == pytest.approx(want)
