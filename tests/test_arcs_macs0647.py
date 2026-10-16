"""Tests for scripts/arcs_macs0647.py: the frozen outcome rule (synthetic summaries; no network)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

from astropy.table import Table

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("arcs_macs0647", _DIR / "arcs_macs0647.py")
am = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(am)


def _run(model, n, aligned, p):
    return {
        "model": model,
        "strong_shear_background": {"n": n, "aligned": aligned, "anti": 0},
        "anti_expected_shape_noise_background": {"p_ge_observed": p},
    }


def test_control_threshold_inclusive_and_empty_fails():
    assert am.control_passes({"n": 10, "aligned": 7})
    assert not am.control_passes({"n": 10, "aligned": 6})
    assert not am.control_passes({"n": 0, "aligned": 0})


def test_primary_decides_when_its_control_passes():
    r = am.outcome(_run("p", 40, 30, 0.01), _run("s", 40, 30, 0.9))
    assert r == {"deciding_model": "p", "outcome": "replicates", "p_ge_observed": 0.01}
    r = am.outcome(_run("p", 40, 30, 0.05), _run("s", 40, 30, 0.001))
    assert r["outcome"] == "does not replicate" and r["deciding_model"] == "p"


def test_second_decides_only_after_primary_control_fails():
    r = am.outcome(_run("p", 40, 20, 0.001), _run("s", 40, 29, 0.2))
    assert r["deciding_model"] == "s" and r["outcome"] == "does not replicate"
    r = am.outcome(_run("p", 40, 20, 0.001), _run("s", 40, 20, 0.001))
    assert r == {"deciding_model": None, "outcome": "no result", "p_ge_observed": None}


def test_match_class_by_label_and_position():
    rows = Table({"label": [1, 2], "ra": [10.0, 10.01], "dec": [70.0, 70.0]})
    other = Table(
        {
            "label": [2, 3],
            "ra": [10.0 + 0.2 / 3600 / 0.342, 10.02],
            "dec": [70.0, 70.0],
            "orientation_class": ["aligned", "anti"],
        }
    )
    assert am.match_class(rows, other, None) == ["not selected", "aligned"]
    assert am.match_class(rows, other, 0.3) == ["aligned", "not selected"]
    assert am.match_class(rows, other[:0], 0.3) == ["not selected", "not selected"]


def test_validate_guard():
    assert am.validate_passes({"kappa_map": {"median_rel_diff": 2e-5}})
    assert not am.validate_passes({"kappa_map": {"median_rel_diff": 0.02}})
