"""Offline tests for scripts/m0717_crossmodel.py (frozen rules; no network)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import m0717_crossmodel as mx  # noqa: E402

SUMMARY = ROOT / "results/lens_consistency/macs0717_crossmodel/summary.json"


def test_combine():
    assert mx.combine([]) == ""
    assert mx.combine(["model-dependent (CATS only)", "both models"]) == "both models"
    assert (
        mx.combine(["model-dependent (CATS only)", "catalogued image misplaced by CATS"])
        == "catalogued image misplaced by CATS"
    )
    assert mx.combine(["model-dependent (CATS only)"]) == "model-dependent (CATS only)"


def test_outcome():
    assert mx.outcome("recovered", "both models") == "no_flag"
    assert mx.outcome("no_reference", "") == "no_flag"
    assert mx.outcome("absent", "") == "untestable"
    assert mx.outcome("absent", "both models") == "lensing_residual_input"
    assert mx.outcome("confused", "both models") == "second_model_supported"
    assert mx.outcome("inconsistent_reference", "both models") == "second_model_supported"
    assert mx.outcome("absent", "model-dependent (CATS only)") == "explained_model_dependent"
    assert mx.outcome("confused", "catalogued image misplaced by CATS") == (
        "explained_model_dependent"
    )


def test_final():
    # independent Step 6 record (68 on a critical curve), within 1"
    assert mx.final("lensing_residual_input", "absent", "68", 109.395425, 37.739381) == (
        "explained_step6"
    )
    # Step 6 CATS-only row now supported by a second model: no independent explanation
    assert mx.final("second_model_supported", "confused", "45", 109.384136, 37.743892) == (
        "to_vet_candidate"
    )
    # a confused copy-system row without its own record
    assert mx.final("second_model_supported", "confused", "12", 109.0, 37.0) == "explained_step6"
    assert mx.final("lensing_residual_input", "absent", "99", 109.0, 37.0) == "to_vet_candidate"
    # system 65: only its inconsistent_reference copies are covered
    assert mx.final("lensing_residual_input", "absent", "65", 109.0, 37.0) == "to_vet_candidate"
    assert mx.final("second_model_supported", "inconsistent_reference", "65", 109.0, 37.0) == (
        "explained_step6"
    )
    assert mx.final("explained_model_dependent", "absent", "3", 109.0, 37.0) == (
        "explained_model_dependent"
    )
    # the record lookup needs the same system and a position within 1"
    assert mx.step6_record("68", 109.395425 + 2.0 / 3600, 37.739381) is None
    assert mx.step6_record("67", 109.395425, 37.739381) is None


def test_tracked_summary_consistent():
    s = json.loads(SUMMARY.read_text())
    used = s["models_used"]
    for name, g in s["gates"].items():
        assert g["passes"] == mx.mc.gate(g["reproduced"], g["tested"])
        assert (name in used) == g["passes"]
    counts: dict[str, int] = {}
    for r in s["rows"]:
        comb = mx.combine([r["models"][n]["d041"] for n in used if n in r["models"]])
        assert r["combined"] == comb
        assert r["outcome"] == mx.outcome(r["forced_class"], comb)
        assert r["final"] == mx.final(r["outcome"], r["forced_class"], r["system"], *r["radec"])
        counts[r["final"]] = counts.get(r["final"], 0) + 1
    assert s["outcome"] == (counts if used else "tool not validated")
    assert s["n_rows"] == len(s["rows"])


def test_tracked_post_hoc_consistent():
    s = json.loads(SUMMARY.read_text())
    p = s["post_hoc_arcs_z"]
    assert p["models_used"] == [n for n, g in s["gates"].items() if g["secondary_arcs_z"]["passes"]]
    frozen = {(r["system"], *r["radec"]): r for r in s["rows"]}
    counts: dict[str, int] = {}
    for r in p["rows"]:
        models = frozen[(r["system"], *r["radec"])]["models"]
        comb = mx.combine([models[n]["d041"] for n in p["models_used"] if n in models])
        assert r["combined"] == comb
        assert r["outcome"] == mx.outcome(r["forced_class"], r["combined"])
        assert r["final"] == mx.final(r["outcome"], r["forced_class"], r["system"], *r["radec"])
        counts[r["final"]] = counts.get(r["final"], 0) + 1
    assert p["outcome"] == counts
