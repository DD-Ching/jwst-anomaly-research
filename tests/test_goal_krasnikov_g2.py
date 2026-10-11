"""Offline tests for scripts/goal_krasnikov_g2.py (gap G2: ANEC on the Krasnikov tube's rays)."""

import importlib.util
import sys
from dataclasses import replace
from pathlib import Path

import pytest

_DIR = Path(__file__).resolve().parents[1] / "scripts"
_spec = importlib.util.spec_from_file_location("goal_krasnikov_g2", _DIR / "goal_krasnikov_g2.py")
g2 = importlib.util.module_from_spec(_spec)
sys.modules["goal_krasnikov_g2"] = g2
_spec.loader.exec_module(g2)

SMALL = replace(g2.P, d_len=60.0, rho_max=10.0, n_rays=8)


def test_outgoing_lines_are_shear_free_null_pregeodesics_with_zero_null_ricci():
    out = g2.outgoing_lines()
    assert out["g_xx"] == "0"
    assert out["pregeodesic"]
    assert out["R_xx"] == "0"  # for an arbitrary k(u, x, rho): front and end caps included
    assert set(out["screen_tidal_R_AxBx"].values()) == {"0"}


def test_region_with_k_of_x_only_is_flat():
    assert g2.flat_where_k_depends_on_x_only()["nonzero_riemann_components"] == []


def test_outgoing_speed_never_exceeds_light():
    for k in (1.0, 0.3, 1e-3):
        assert g2.max_outgoing_speed(k, n=4000) <= 1.0


def test_compact_step_is_exact_outside_its_transition():
    assert g2.theta(-0.1, 1.0) == 0.0 and g2.theta(1.0, 1.0) == 1.0 and g2.theta(5.0, 1.0) == 1.0
    assert 0.0 < g2.theta(0.5, 1.0) < 1.0


def test_core_rays_after_completion_stay_flat_and_early_rays_do_not():
    late = g2.incoming_rays(SMALL)
    assert late["max_dk_drho_or_du_along_rays"] == 0.0
    assert late["transit_time_max_eps"] < 0.5 * late["exterior_light_time_eps"]  # the shortcut
    early = g2.early_rays(SMALL)  # positive control: launched before completion, meets the front
    assert early["max_dk_drho_or_du_along_rays"] > 0.1


def test_axis_regularity_of_er_profiles():
    # owner's value: ER tanh radial factor, eta = 0.999, R = 10 eps -> k_rho(0) = 7.87544e-33 / eps
    assert float(g2.radial_derivative_at_axis_tanh(0.999, 10.0)) == pytest.approx(
        7.87544e-33, rel=1e-5
    )
    assert abs(g2.front_axis_derivative(0.5)) > 0.1  # front crossing the axis: kink
    assert g2.front_axis_derivative(2.0) == 0.0  # after the front: smooth


def test_axial_normalisation_and_anec_formula():
    ax = g2.axis_null_ricci()
    assert ax["R_ll_on_axis"] == "-(k0(x) + 1)*k2(x)/k0(x)**2"  # for k smooth at the axis
    n = g2.axial_normalisation()
    assert n["g_ll"] == "0" and n["killing_energy"] == "1"
    assert n["R_ll_affine"] == "-4*k2/(k0 + 1)"
    assert n["anec_integrand_per_dx"] == "-2*k2"


def test_smooth_axis_examples():
    ex = g2.smooth_axis_examples(d_len=60.0)
    a, b = ex["exp_rho2"], ex["exp_rho4"]
    assert a["anec_minus_2_int_k2"] < 0
    assert a["anec_direct_rho_1e-5R"] == pytest.approx(a["anec_minus_2_int_k2"], rel=1e-4)
    assert b["anec_minus_2_int_k2"] == 0.0 and abs(b["anec_direct_rho_1e-5R"]) < 1e-6
    assert b["R_ll_unit_at_rho_R_over_2_mid_tube"] < 0  # k2 = 0 does not mean a flat core


def test_w_never_decreases_along_causal_vectors():
    for step in (g2.theta, g2.theta_tanh):
        assert g2.w_monotone_check(step, n=1500, p=SMALL)["min_normalised_dw"] >= 0.0
