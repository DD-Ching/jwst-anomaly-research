"""Tests for scripts/arcs_anti_ledger.py: class selection, row dropping, frozen outcome."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest
from astropy.table import Table

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("arcs_anti_ledger", _DIR / "arcs_anti_ledger.py")
al = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(al)


def _table():
    return Table(
        {
            "field": ["F", "F", "F", "F"],
            "label": ["1", "2", "3", "4"],
            "gate_passed": [True, True, True, False],
            "deblend_class": ["single", "single", "two-knot", "single"],
            "screen_class": ["anti", "aligned", "anti", "anti"],
            "p_anti_row": [0.1, 0.2, 0.5, 0.5],
        }
    )


def test_class_keeps_gate_passed_non_two_knot():
    assert list(al.class_rows(_table())["label"]) == ["1", "2"]


def test_variant_drops_row_and_probability():
    c = al.class_rows(_table())
    full = al.variant(c)
    assert (full["n"], full["anti"], full["expected"]) == (2, 1, 0.3)
    assert full["p_ge_k"] == pytest.approx(1 - 0.9 * 0.8)
    d = al.variant(c, ("F 1",))
    assert (d["n"], d["anti"], d["p_ge_k"]) == (1, 0, 1.0)
    o = al.variant(c, ("F 1",), drop_probs=False)
    assert (o["n"], o["anti"], o["expected"]) == (2, 0, 0.3)
    with pytest.raises(ValueError):
        al.variant(c, ("F 9",))


def test_decide_threshold():
    assert al.decide(0.05).startswith("single-component excess not supported")
    assert al.decide(0.0499).startswith("single-component excess remains")
    assert al.decide(float("nan")) == "undetermined"


def test_stored_census_reproduces_and_runs():
    res = al.run(Table.read(al.CENSUS))
    assert res["reproduce"]["anti"] == 3
    assert res["primary"]["anti"] == 1
