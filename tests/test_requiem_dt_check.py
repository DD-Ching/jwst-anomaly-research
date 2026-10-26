"""Rules of scripts/requiem_dt_check.py (plan results/lens_consistency/requiem_dt/PLAN.md)."""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import requiem_dt_check as r  # noqa: E402


def test_observed_mapping_signs():
    d = r.observed_draws(n=20_000)
    ba, ca = r.summarize(d["ba"]), r.summarize(d["ca"])
    assert math.isclose(ba[0], -114.0, abs_tol=1.5)
    assert math.isclose(ca[0], -131.0, abs_tol=2.5)


def test_delay_z_zero_for_equal_values():
    m = {"ba": (-50.0, 5.0, 5.0)}
    assert r.delay_z(m, "ba", (-50.0, 10.0, 10.0)) == 0.0
    assert math.isclose(r.delay_z(m, "ba", (-80.0, 10.0, 10.0)), -30 / math.hypot(10, 5))


def test_ratio_z_zero_when_model_matches():
    robs = {"median": 1.5, "sigma_ln": 0.2}
    m = {"ba": (-40.0, 4.0, 4.0), "ca": (-60.0, 6.0, 6.0)}
    assert math.isclose(r.ratio_z(m, robs), 0.0, abs_tol=1e-12)


def test_outcome_rules():
    rows = {"a": {"ratio": 4.0, "ba": 0.0, "ca": 0.0}, "b": {"ratio": -3.5, "ba": 4.0, "ca": 0.0}}
    assert r.outcome(rows, True) == (["ratio"], "F")
    assert r.outcome(rows, False) == ([], "U")
    rows["b"]["ratio"] = 2.0
    assert r.outcome(rows, True) == ([], "N")
