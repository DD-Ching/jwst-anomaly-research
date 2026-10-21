"""Offline tests for scripts/m1149_crossmodel.py (frozen outcome and gate; no network)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import m1149_crossmodel as mc  # noqa: E402

SUMMARY = ROOT / "results/lens_consistency/macs1149_crossmodel/summary.json"


def test_outcome_rule():
    assert mc.outcome(False, "absent", "both models") == "untestable"
    assert mc.outcome(True, "absent", "both models") == "lensing_residual_input"
    assert mc.outcome(True, "absent", "model-dependent (CATS only)") == "explained_model_dependent"
    assert (
        mc.outcome(True, "absent", "catalogued image misplaced by CATS")
        == "explained_model_dependent"
    )
    assert mc.outcome(True, "recovered", "both models") == "no_flag"
    assert mc.outcome(True, "no_reference", "both models") == "no_flag"


def test_gate():
    assert mc.gate(8, 10)
    assert not mc.gate(7, 10)
    assert not mc.gate(0, 0)


def test_tracked_summary_consistent():
    s = json.loads(SUMMARY.read_text())
    assert s["gate"]["passes"] == mc.gate(s["gate"]["reproduced"], s["gate"]["tested"])
    counts: dict[str, int] = {}
    for r in s["rows"]:
        d041 = r.get("sharon", {}).get("d041", "")
        assert r["outcome"] == mc.outcome(r["inside_footprint"], r["forced_class"], d041)
        counts[r["outcome"]] = counts.get(r["outcome"], 0) + 1
    assert s["outcome"] == counts
    assert len(s["rows"]) == 16
