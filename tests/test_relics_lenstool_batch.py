"""Frozen rules of scripts/relics_lenstool_batch.py (H_K plan), offline."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import relics_lenstool_batch as rb  # noqa: E402


def test_classes_follow_plan_thresholds():
    assert rb.classify(1.25) == "offset"
    assert rb.classify(1.50) == "offset"
    assert rb.classify(1.00) == "near_1"
    assert rb.classify(0.85) == "near_1"
    assert rb.classify(1.20) == "intermediate"
    assert rb.classify(0.80) == "intermediate"


def test_ratio_stats_uses_reference_region_and_ignores_nan():
    ref = np.array([10.0, 9.0, 6.0, 1.0, np.nan])
    lt = 1.5 * ref
    lt[3] = 100.0  # outside the >= 0.5 max region: must not enter R
    s = rb.ratio_stats(lt, ref)
    assert np.isclose(s["R"], 1.5)
    assert s["region_arcsec2"] == 3.0
    assert s["overlap_arcsec2"] == 4.0
    assert s["region_fraction"] == 0.75


def _rows(classes_r):
    return {c: {"R": r, "class": rb.classify(r)} for c, r in classes_r.items()}


FACTORS = {"macs0553": 0.665, "a3192": 0.692, "rxc2211": 1.017, "cl0152": 1.045}


def test_gate_v_and_h_k_outcomes():
    good = {"macs0553": 1.5, "a3192": 1.45, "rxc2211": 0.98, "cl0152": 0.96}
    rows = _rows(good | {"rxc0600": 1.4, "macs0417": 1.5, "whl0137": 1.0, "spt0615": 1.05})
    gate = rb.gate_v(rows, FACTORS)
    assert gate["passed"]
    assert rb.h_k(rows, gate)["outcome"] == "supported"
    rows["spt0615"] = {"R": 1.4, "class": "offset"}
    assert rb.h_k(rows, gate)["outcome"] == "rejected"
    rows["spt0615"] = {"R": 1.2, "class": "intermediate"}
    assert rb.h_k(rows, gate)["outcome"].startswith("undetermined")
    bad = _rows(good | {"rxc2211": 1.3})
    gate_bad = rb.gate_v(bad, FACTORS)
    assert not gate_bad["passed"]
    assert rb.h_k(bad, gate_bad)["outcome"].startswith("undetermined")


def test_gate_v_tolerance_applies_within_class():
    rows = _rows({"macs0553": 1.26, "a3192": 1.45, "rxc2211": 0.98, "cl0152": 0.96})
    gate = rb.gate_v(rows, FACTORS)  # 1/1.26 = 0.794, 0.129 from 0.665
    assert gate["rows"]["macs0553"]["class_ok"] and not gate["rows"]["macs0553"]["tol_ok"]
    assert not gate["passed"]
