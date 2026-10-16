"""Tests for scripts/arcs_macs1423.py: the frozen outcome rule (synthetic summaries; no network)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("arcs_macs1423", _DIR / "arcs_macs1423.py")
am = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(am)


def _run(n, aligned, p):
    return {
        "strong_shear_background": {"n": n, "aligned": aligned, "anti": 0},
        "anti_expected_shape_noise_background": {"p_ge_observed": p},
    }


def test_control_failure_gives_no_result_whatever_p():
    assert am.outcome(_run(40, 27, 0.001)) == {"outcome": "no result", "p_ge_observed": None}
    assert am.outcome(_run(0, 0, 0.001))["outcome"] == "no result"


def test_replication_threshold_is_strict():
    assert am.outcome(_run(40, 28, 0.049)) == {"outcome": "replicates", "p_ge_observed": 0.049}
    assert am.outcome(_run(40, 28, 0.05))["outcome"] == "does not replicate"
    assert am.outcome(_run(10, 7, 0.9))["outcome"] == "does not replicate"
