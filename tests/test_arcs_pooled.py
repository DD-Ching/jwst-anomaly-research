"""Tests for scripts/arcs_pooled.py: frozen outcome rule, pooling, post hoc S1063 comparison."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, _DIR / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ap = _load("arcs_pooled")


@pytest.mark.parametrize(
    ("pa", "pb", "want"),
    [
        (0.01, 0.04, "pooled excess"),
        (0.05, 0.2, "no pooled excess"),
        (0.01, 0.05, "mixed"),
        (0.2, 0.01, "mixed"),
    ],
)
def test_outcome_rule(pa, pb, want):
    assert ap.outcome(pa, pb) == want


def test_pooled_sums_and_poisson():
    stats = {"a": {"anti": 2, "expected": 0.5}, "b": {"anti": 0, "expected": 0.5}}
    r = ap.pooled(stats, ("a", "b"))
    assert r["anti"] == 2 and r["expected"] == 1.0
    assert r["p_ge_observed"] == pytest.approx(1 - 2 / 2.718281828459045, rel=1e-6)


def test_variant_b_drops_only_macs1423():
    assert set(ap.VARIANTS["A"]) - set(ap.VARIANTS["B"]) == {"MACS1423"}


def test_s1063_vs_pooled_scales_and_conditions():
    r = ap.s1063_vs_pooled({"anti": 1, "expected": 1.0}, {"anti": 2, "expected": 1.0})
    assert r["pool_ratio"] == 2.0 and r["scaled_expected"] == 2.0
    assert r["p_conditional"] == pytest.approx(1 - 0.5**3)  # P(>= 1 of 3 at p = 0.5)
