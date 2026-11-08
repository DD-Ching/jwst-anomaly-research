"""Tests for scripts/arcs_pooled_census.py: frozen outcome rule, row probabilities, rows."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, _DIR / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


pc = _load("arcs_pooled_census")


def test_decide_gates_then_outcome():
    assert pc.decide(11, 20, 0.5).startswith("undetermined")
    assert pc.decide(12, 30, 0.5).startswith("undetermined")  # pass fraction 0.4 < 0.5
    assert pc.decide(0, 0, 0.5).startswith("undetermined")
    assert pc.decide(60, 104, 0.05) == "anti excess explained by clump orientation (exploratory)"
    assert pc.decide(60, 104, 0.049) == "excess remains under the clump-aware null"
    assert pc.decide(60, 104, float("nan")).startswith("undetermined")


def test_row_probs_follow_273_rules():
    p_cp, p_sn, p = pc.row_probs(0.5, "two-knot", True, 0.3)
    assert p == p_cp and p_cp > p_sn
    assert pc.row_probs(0.5, "single", True, 0.3)[2] == p_sn
    mix = pc.row_probs(0.5, "single", False, 0.3)[2]
    assert mix == pytest.approx(0.3 * p_cp + 0.7 * p_sn)
    # g_min is rounded to 2 decimals as arcs_pooled.field_stats
    assert pc.row_probs(0.5004, "single", True, 0.0)[1] == p_sn


def test_rows_reproduce_pooled_tally():
    rows = pc.load_rows()
    assert len(rows) == 104
    assert int(np.sum(rows["orientation_class"] == "anti")) == 6
    lam = sum(pc.row_probs(g, "single", True, 0.0)[1] for g in rows["g_min"])
    assert round(lam, 3) == pc.POOLED_B_EXPECTED
    assert {
        (str(r["field"]), str(r["label"])) for r in rows if r["orientation_class"] == "anti"
    } >= pc.PERSIST


def test_image_urls():
    assert pc.image_url("WHL0137").endswith(
        "jw02282/L3/t/o010/jw02282-o010_t001_nircam_clear-f200w_i2d.fits"
    )
    assert set(pc.IMAGES) == set(pc.FIELDS_B)


def test_posthoc_class_split():
    from astropy.table import Table

    t = Table(
        {
            "gate_passed": [True, True, True, False],
            "deblend_class": ["two-knot", "single", "single", "undecided"],
            "screen_class": ["anti", "anti", "aligned", "aligned"],
            "p_anti_row": [0.2, 0.1, 0.1, 0.3],
        }
    )
    ph = pc.posthoc(t)
    assert ph["two_knot_gate_passed"] == {
        "n": 1,
        "anti": 1,
        "expected": 0.2,
        "p_ge_k": pytest.approx(0.2),
        "p_le_k": pytest.approx(1.0),
    }
    assert ph["other_gate_passed"]["anti"] == 1 and ph["other_gate_passed"]["n"] == 2
    assert ph["gate_failed"]["anti"] == 0 and ph["gate_failed"]["p_ge_k"] == pytest.approx(1.0)


def test_persisting_tail_drops_non_persisting_anti():
    keys = [("WHL0137", "2194"), ("MACS0417", "1993"), ("X", "1")]
    probs = [0.5, 0.9, 0.5]
    anti = keys[:2]
    # 1993 leaves count and null: P(>= 1) over two rows of 0.5
    assert pc.persisting_tail(keys, probs, anti) == pytest.approx(0.75)


def test_finite():
    assert pc._finite(float("nan")) is None and pc._finite(0.3) == 0.3
