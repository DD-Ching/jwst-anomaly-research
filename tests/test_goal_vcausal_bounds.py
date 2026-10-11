"""Offline tests for scripts/goal_vcausal_bounds.py (positive control: published bounds)."""

from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import pytest

_DIR = Path(__file__).resolve().parents[1] / "scripts"
_spec = importlib.util.spec_from_file_location(
    "goal_vcausal_bounds", _DIR / "goal_vcausal_bounds.py"
)
vb = importlib.util.module_from_spec(_spec)
sys.modules["goal_vcausal_bounds"] = vb
_spec.loader.exec_module(vb)


def test_formula_reproduces_the_three_published_bounds():
    for row in vb.calibration():
        assert row["ratio"] == pytest.approx(1.0, abs=0.02), row["experiment"]


def test_bound_limits():
    assert vb.speed_bound(1e-3, 0.0, 0.0) == pytest.approx(math.sqrt(1 + (1 - 1e-6) / 1e-6))
    assert vb.speed_bound(0.0, 0.0, 0.0) == math.inf
    # a larger offset or residual can only weaken the bound
    assert vb.speed_bound(1e-6, 1e-3, 1e-6) > vb.speed_bound(1e-5, 1e-3, 1e-6)
    with pytest.raises(ValueError):
        vb.speed_bound(1.0, 0.0, 0.0)


def test_nist_2015_cannot_compete_with_dedicated_tests():
    best = max(r["best_bound"] for r in vb.nist_best_case())
    assert best < 100  # vs ~5e6 c already published (Cocciaro et al. 2018)


def test_requirements_close_at_target():
    v = vb.required_speed(20, 1.0)
    assert v == pytest.approx(6.31e8, rel=1e-3)
    req = vb.requirements(v)
    assert req["check_bound_at_budget"] == pytest.approx(v, rel=1e-6)
