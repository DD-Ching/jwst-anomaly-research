"""Rules of scripts/refsdal_sx_check.py (plan results/lens_consistency/refsdal_sx/PLAN.md)."""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import refsdal_sx_check as r  # noqa: E402


def test_observed_sigmas():
    s_dt, s_ln = r.obs_sigma()
    assert math.isclose(s_dt, 5.575)
    assert math.isclose(s_ln, 0.5 * math.log(0.354 / 0.272))


def test_z_zero_at_observed_values():
    assert r.delay_z((r.OBS_DT[0], 10.0, 10.0)) == 0.0
    assert math.isclose(r.ratio_z((r.OBS_MU[0], 0.05, 0.05)), 0.0, abs_tol=1e-12)
    s_dt, _ = r.obs_sigma()
    assert math.isclose(r.delay_z((r.OBS_DT[0] - 30, 4.0, 2.0)), 30 / math.hypot(s_dt, 3.0))


def test_outcome_rules():
    ok = {"flagged": False, "both_given": True}
    bad = {"flagged": True, "both_given": True}
    dt_only = {"flagged": True, "both_given": False}
    assert r.outcome({"a": bad, "b": ok}) == "N"
    assert r.outcome({"a": bad, "b": bad, "c": {"flagged": False, "both_given": False}}) == "F"
    assert r.outcome({"a": dt_only}) == "U"


def test_flag_either_quantity():
    m = {"dt": (r.OBS_DT[0], 1.0, 1.0), "mu": (r.OBS_MU[0] * 3, 0.01, 0.01)}
    assert r.row(m)["flagged"]
