"""Offline tests for scripts/goal_khronometric.py (goal round 9, P3)."""

import importlib.util
import sys
from pathlib import Path

import sympy as sp

_DIR = Path(__file__).resolve().parents[1] / "scripts"
_spec = importlib.util.spec_from_file_location("goal_khronometric", _DIR / "goal_khronometric.py")
kh = importlib.util.module_from_spec(_spec)
sys.modules["goal_khronometric"] = kh
_spec.loader.exec_module(kh)


def _zero(expr: str) -> bool:
    return sp.simplify(sp.sympify(expr)) == 0


def test_blas_sibiryakov_controls():
    c = kh.bs_controls()
    assert _zero(c["partial_fraction_residual"])
    assert _zero(c["instantaneous_trace_position_space"])
    assert all(_zero(d) for d in c["instantaneous_divergence_position_space"])
    # gradient source: no 1/p^2 left, so no long-range instantaneous term
    assert "/(px**2" not in c["gradient_source_exchange"].replace(" ", "")
    assert "Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]])" == c["ea_instantaneous_sum"]


def test_khronometric_matter_response_equals_ea_for_c13_zero():
    r = kh.khronometric_vs_ea()
    assert r["n_solutions"] == 1
    assert all(_zero(v) for v in r["difference_c13_0"].values())
    for comp in ("xx", "yy", "zz", "xy"):
        assert _zero(r["difference_general"][comp])


def test_c13_nonzero_instantaneous_part():
    r = kh.c13_check()
    assert all(_zero(v) for v in r["frame_residuals_EA"].values())
    assert all(_zero(v) for v in r["frame_residuals_khronometric"].values())
    assert all(_zero(r["coefficient_difference_kh_minus_EA"][n]) for n in "abcdt")
    assert _zero(r["vector_difference_minus_owner_form"])
    assert all(_zero(v) for v in r["khronometric_d_dc1_at_fixed_c13_c14"].values())
    assert all(_zero(v) for v in r["split_residual"])
    bare = r["bare_k2_factors"]
    for name in ("EA GR limit", "EA c13=0", "khronometric c13=0", "EA c13=3/20", "EA c13=1/1000"):
        assert bare[name] == [], name
    for c13 in ("3/20", "1/1000"):
        assert bare[f"khronometric c13={c13}"] == ["kx**2 + ky**2 + kz**2"]
        assert bare[f"khronometric - I c13={c13}"] == []
    assert all(_zero(v["static_khronometric_minus_EA"]) for v in r["pole_bookkeeping"].values())


def test_position_space_instantaneous_kernel():
    r = kh.position_space_check(n=96, box=12.0)
    assert all(p["max_rel_error"] < 0.1 for p in r["points"])


def test_link_budget_c13_scaling():
    r = kh.link_budget_c13()
    assert r["angular_factor"] == 36.0
    rows = {(x["distance_ly"], x["receive_time_s"]): x for x in r["rows"]}
    # kappa needed scales as r^5 / sqrt(tau)
    ratio = rows[(100.0, 1.0)]["kappa_required"] / rows[(20.0, 1.0)]["kappa_required"]
    assert abs(ratio / 5**5 - 1) < 1e-9


def test_ppn_identity_and_window_bound():
    r = kh.ppn_identity()
    assert _zero(r["identity_residual"])
    w = kh.ppn_window_check()
    assert abs(w["coef_cS2"] - 6.3e-6) < 1e-12
    for row in w["rows"]:
        assert row["kappa_allowed_PPN"] < row["kappa_required"] * 1e-40
        assert row["kappa_cap_perturbative_only"] < row["kappa_required"] * 1e-20
        assert row["alpha2_v2_phi_floor"] > 1e20  # the needed points lie outside control
