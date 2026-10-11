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
