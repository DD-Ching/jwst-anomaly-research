"""Tests for scripts/arcs_vet_census.py frozen rules (no network)."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np

_DIR = Path(__file__).resolve().parents[1] / "scripts"
_spec = importlib.util.spec_from_file_location("arcs_vet_census", _DIR / "arcs_vet_census.py")
vs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(vs)


def test_gate_separation_clips_and_rounds():
    assert vs.gate_separation(0.05) == 0.10
    assert vs.gate_separation(0.147) == 0.15
    assert vs.gate_separation(0.9) == 0.30


def test_row_p_anti():
    assert vs.row_p_anti("two-knot", True, 0.1, 0.01, 0.5) == 0.1
    assert vs.row_p_anti("single", True, 0.1, 0.01, 0.5) == 0.01
    assert vs.row_p_anti("undecided", True, 0.1, 0.01, 0.5) == 0.01
    assert np.isclose(vs.row_p_anti("two-knot", False, 0.1, 0.01, 0.25), 0.25 * 0.1 + 0.75 * 0.01)


def test_clopper_pearson_upper():
    assert vs.clopper_pearson_upper(5, 5) == 1.0
    assert np.isclose(vs.clopper_pearson_upper(0, 10), 1 - 0.025 ** (1 / 10))


def test_outcome():
    assert vs.outcome(float("nan")).startswith("undetermined")
    assert vs.outcome(0.05).startswith("anti excess explained")
    assert vs.outcome(0.049).startswith("excess remains")
